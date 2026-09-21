#!/usr/bin/env python3
"""Retrieve precomputed microexon code-model predictions for variants in a VCF."""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, TextIO

try:
    import pysam
except ImportError as exc:  # pragma: no cover - exercised by users without dependencies
    raise SystemExit(
        "pysam is required. Install it with: "
        "conda install -c conda-forge -c bioconda pysam=0.24.1"
    ) from exc


OUTPUT_COLUMNS = [
    "variant",
    "geneName",
    "var_type",
    "event",
    "ME-ID",
    "MIC_group",
    "MIC_dist",
    "var_MIC_region",
    "mt_code_score",
    "wt_code_score",
    "logit_scaling_factor",
    "delta_code_model",
    "assembly",
    "VCF_ID",
    "prediction_status",
    "ref_ok",
]

EVENT_COLUMNS = [
    "event",
    "chrom",
    "strand",
    "upIntStart",
    "upIntEnd",
    "dnIntStart",
    "dnIntEnd",
    "geneName",
    "lengthDiff",
    "group",
]

INDEL_COLUMNS = [
    "variant",
    "event",
    "var_MIC_region",
    "MIC_dist",
    "mt_code_score",
    "wt_code_score",
    "logit_scaling_factor",
    "delta_code_model",
    "logit_scaling_factor_all_vars",
    "delta_code_model_all_vars",
]

SCALING_COLUMNS = [
    "event",
    "logit_scaling_factor",
    "logit_scaling_factor_all_vars",
]

DELTA_SCALING_COLUMNS = {
    "all-variants": (
        "logit_scaling_factor_all_vars",
        "delta_code_model_all_vars",
    ),
    "main-analysis": ("logit_scaling_factor", "delta_code_model"),
}

ISM_COLUMNS = [
    "event",
    "region",
    "rel_pos",
    "ref",
    "alt",
    "mt_score",
    "delta_logit_score",
    "transformed_delta_scores",
]

DNA = set("ACGTN")
BIN_SIZE = 1_000


def normalize_chrom(chrom: str) -> str:
    """Normalize common GRCh38 contig names to UCSC-style names."""
    chrom = str(chrom).strip()
    if chrom.lower().startswith("chr"):
        chrom = chrom[3:]
    if chrom.upper() in {"M", "MT"}:
        return "chrM"
    return f"chr{chrom}"


def canonical_group(group: str) -> str:
    return "new_mic" if group == "novel_mic" else group


def is_simple_snv(ref: str, alt: str) -> bool:
    return len(ref) == 1 and len(alt) == 1 and ref in DNA - {"N"} and alt in DNA - {"N"}


def is_sequence_allele(allele: str) -> bool:
    return bool(allele) and set(allele) <= DNA


def classify_variant(ref: str, alt: str) -> str:
    if is_simple_snv(ref, alt):
        return "SNV"
    if is_sequence_allele(ref) and is_sequence_allele(alt) and len(ref) != len(alt):
        return "Indel"
    return "Unsupported"


def logit2(probability: float) -> float:
    probability = float(probability)
    epsilon = 1e-12
    probability = min(max(probability, epsilon), 1 - epsilon)
    return math.log2(probability / (1 - probability))


def require_columns(actual: list[str] | None, expected: list[str], source: Path) -> None:
    missing = [column for column in expected if column not in (actual or [])]
    if missing:
        raise ValueError(f"{source} is missing required columns: {', '.join(missing)}")


@dataclass(frozen=True)
class Event:
    event: str
    chrom: str
    strand: str
    length_diff: int
    c1: int
    c2: int
    me_start: int
    me_end: int
    gene_name: str
    microexon_group: str
    microexon_id: str


@dataclass(frozen=True)
class Interval:
    event: str
    region: str
    start0: int
    end0: int


@dataclass
class QueryHit:
    chrom: str
    pos: int
    record_id: str
    ref: str
    alt: str
    variant: str
    event: str
    region: str
    relative_position: int
    mic_distance: int | str
    variant_type: str
    prediction: dict[str, str] | None = None
    prediction_status: str = "NA_NOT_PRECOMPUTED"
    reference_check: str = "NOT_REQUESTED"


class PredictionData:
    def __init__(self, data_dir: Path, delta_scaling: str = "all-variants"):
        if delta_scaling not in DELTA_SCALING_COLUMNS:
            choices = ", ".join(DELTA_SCALING_COLUMNS)
            raise ValueError(f"Unknown delta scaling mode {delta_scaling!r}; choose {choices}")
        self.data_dir = Path(data_dir)
        self.delta_scaling = delta_scaling
        (
            self.scaling_factor_column,
            self.delta_code_model_column,
        ) = DELTA_SCALING_COLUMNS[delta_scaling]
        self.ism_path = self.data_dir / "MIC_ISM.gz"
        self.wt_path = self.data_dir / "MIC_ISM_wt.gz"
        self.events_path = self.data_dir / "microexon_events.hg38.csv.gz"
        self.indels_path = self.data_dir / "MIC_indels.csv.gz"
        self.scaling_path = self.data_dir / "MIC_event_scaling.csv.gz"

        for path in (
            self.ism_path,
            self.wt_path,
            self.events_path,
            self.indels_path,
            self.scaling_path,
        ):
            if not path.is_file():
                raise FileNotFoundError(f"Required data file not found: {path}")

        self.wt_by_event = self._load_wild_type_scores()
        self.events, self.interval_bins = self._load_events()
        self.indels_by_variant = self._load_indels()
        self.scale_by_event = self._load_scaling()
        self._validate_coverage()

    def _load_wild_type_scores(self) -> dict[str, str]:
        scores: dict[str, str] = {}
        with gzip.open(self.wt_path, "rt", newline="") as handle:
            reader = csv.DictReader(handle)
            require_columns(reader.fieldnames, ["event", "wt_score"], self.wt_path)
            for row in reader:
                event = row["event"]
                score = row["wt_score"]
                if event in scores and scores[event] != score:
                    raise ValueError(f"Conflicting wild-type scores for {event}")
                scores[event] = score
        return scores

    def _load_events(self) -> tuple[dict[str, Event], dict[str, dict[int, list[Interval]]]]:
        events: dict[str, Event] = {}
        interval_bins: dict[str, dict[int, list[Interval]]] = defaultdict(lambda: defaultdict(list))

        with gzip.open(self.events_path, "rt", newline="") as handle:
            reader = csv.DictReader(handle)
            require_columns(reader.fieldnames, EVENT_COLUMNS, self.events_path)
            for row in reader:
                event_id = row["event"]
                chrom = normalize_chrom(row["chrom"])
                strand = row["strand"]
                length_diff = int(row["lengthDiff"])

                c1 = int(row["upIntStart"])
                c2 = int(row["dnIntEnd"])
                me_start = int(row["upIntEnd"])
                me_end = int(row["dnIntStart"])
                if strand == "-":
                    c1 = int(row["upIntEnd"])
                    c2 = int(row["dnIntStart"])
                    me_start = int(row["dnIntEnd"])
                    me_end = int(row["upIntStart"])
                elif strand != "+":
                    raise ValueError(f"Unexpected strand for {event_id}: {strand}")

                me_end -= 1
                event = Event(
                    event=event_id,
                    chrom=chrom,
                    strand=strand,
                    length_diff=length_diff,
                    c1=c1,
                    c2=c2,
                    me_start=me_start,
                    me_end=me_end,
                    gene_name=row["geneName"],
                    microexon_group=canonical_group(row["group"]),
                    microexon_id=f"{chrom}_{strand}_{me_start}_{me_end + 1}",
                )
                if event_id in events:
                    raise ValueError(f"Duplicate event metadata row: {event_id}")
                events[event_id] = event

                if strand == "+":
                    windows = {
                        "up": (me_start - 300, me_start - 1),
                        "ex": (me_start, me_start + length_diff - 1),
                        "dn": (me_end + 1, me_end + 300),
                        "c1": (c1 - 300, c1 + 299),
                        "c2": (c2 - 300, c2 + 299),
                    }
                else:
                    windows = {
                        "up": (me_end + 1, me_end + 300),
                        "ex": (me_end - length_diff + 1, me_end),
                        "dn": (me_start - 300, me_start - 1),
                        "c1": (c1 - 300, c1 + 299),
                        "c2": (c2 - 300, c2 + 299),
                    }

                for region, (start0, end0) in windows.items():
                    if end0 < 0:
                        continue
                    interval = Interval(event_id, region, max(0, start0), end0)
                    for bin_id in range(interval.start0 // BIN_SIZE, interval.end0 // BIN_SIZE + 1):
                        interval_bins[chrom][bin_id].append(interval)

        return events, interval_bins

    def _load_indels(self) -> dict[str, list[dict[str, str]]]:
        indels: dict[str, list[dict[str, str]]] = defaultdict(list)
        with gzip.open(self.indels_path, "rt", newline="") as handle:
            reader = csv.DictReader(handle)
            require_columns(reader.fieldnames, INDEL_COLUMNS, self.indels_path)
            for row in reader:
                indels[row["variant"]].append(row)
        return dict(indels)

    def _load_scaling(self) -> dict[str, float]:
        scales: dict[str, float] = {}
        with gzip.open(self.scaling_path, "rt", newline="") as handle:
            reader = csv.DictReader(handle)
            require_columns(reader.fieldnames, SCALING_COLUMNS, self.scaling_path)
            for row in reader:
                event = row["event"]
                scale = float(row[self.scaling_factor_column])
                if event in scales:
                    raise ValueError(f"Duplicate scaling-factor row for {event}")
                if not math.isfinite(scale) or scale <= 0:
                    raise ValueError(
                        f"Invalid {self.scaling_factor_column} for {event}: {scale}"
                    )
                scales[event] = scale
        return scales

    def _validate_coverage(self) -> None:
        wt_events = set(self.wt_by_event)
        metadata_events = set(self.events)
        scale_events = set(self.scale_by_event)
        if metadata_events != wt_events:
            raise ValueError(
                "Event coverage differs between MIC_ISM_wt.gz and microexon metadata: "
                f"WT={len(wt_events)}, metadata={len(metadata_events)}"
            )
        if scale_events != wt_events:
            raise ValueError(
                "Event coverage differs between MIC_ISM_wt.gz and the event scaling table: "
                f"WT={len(wt_events)}, scaling={len(scale_events)}"
            )

    def find_overlaps(self, chrom: str, start0: int, end0: int) -> list[Interval]:
        chrom = normalize_chrom(chrom)
        found: dict[tuple[str, str], Interval] = {}
        bins = self.interval_bins.get(chrom, {})
        for bin_id in range(start0 // BIN_SIZE, end0 // BIN_SIZE + 1):
            for interval in bins.get(bin_id, []):
                if interval.start0 <= end0 and interval.end0 >= start0:
                    found[(interval.event, interval.region)] = interval
        return sorted(found.values(), key=lambda item: (item.event, item.region))

    def stream_snv_predictions(
        self, wanted: set[tuple[str, str, int, str, str]]
    ) -> dict[tuple[str, str, int, str, str], str]:
        predictions: dict[tuple[str, str, int, str, str], str] = {}
        if not wanted:
            return predictions

        with gzip.open(self.ism_path, "rt", newline="") as handle:
            reader = csv.DictReader(handle)
            require_columns(reader.fieldnames, ISM_COLUMNS, self.ism_path)
            for row in reader:
                key = (
                    row["event"],
                    row["region"].lower(),
                    int(row["rel_pos"]),
                    row["ref"].upper(),
                    row["alt"].upper(),
                )
                if key in wanted:
                    predictions[key] = row["mt_score"]
                    if len(predictions) == len(wanted):
                        break
        return predictions


class ReferenceValidator:
    def __init__(self, fasta_path: Path):
        self.fasta_path = Path(fasta_path)
        index_path = Path(f"{self.fasta_path}.fai")
        if not self.fasta_path.is_file():
            raise FileNotFoundError(f"Reference FASTA not found: {self.fasta_path}")
        if not index_path.is_file():
            raise FileNotFoundError(
                f"Reference FASTA index not found: {index_path}. Create it with samtools faidx."
            )
        self.fasta = pysam.FastaFile(str(self.fasta_path))
        self.contig_by_normalized = {
            normalize_chrom(contig): contig for contig in self.fasta.references
        }

    def check(self, chrom: str, pos1: int, ref: str) -> tuple[str, str | None]:
        normalized = normalize_chrom(chrom)
        contig = self.contig_by_normalized.get(normalized)
        if contig is None:
            return "ERROR", f"contig {chrom} is absent from {self.fasta_path}"
        try:
            observed = self.fasta.fetch(contig, pos1 - 1, pos1 - 1 + len(ref)).upper()
        except (ValueError, OSError) as exc:
            return "ERROR", str(exc)
        if observed != ref:
            return "FAIL", f"expected REF {observed}, received {ref}"
        return "PASS", None

    def close(self) -> None:
        self.fasta.close()


def code_model_relative_position(event: Event, region: str, pos0: int) -> int:
    if event.strand == "+":
        if region == "up":
            return pos0 - event.me_start
        if region == "dn":
            return pos0 - event.me_end
        if region == "ex":
            return pos0 - event.me_start + 1
        if region == "c1":
            return pos0 - event.c1
        if region == "c2":
            return pos0 - event.c2
    else:
        if region == "up":
            return event.me_end - pos0
        if region == "dn":
            return event.me_start - pos0
        if region == "ex":
            return event.me_end - pos0 + 1
        if region == "c1":
            return event.c1 - pos0
        if region == "c2":
            return event.c2 - pos0
    raise ValueError(f"Unknown region: {region}")


def st12_mic_distance(event: Event, region: str, pos0: int) -> int:
    """Return MIC_dist using the convention in Supplementary Table 12."""
    if event.strand == "+":
        if region == "up":
            return pos0 - event.me_start
        if region == "dn":
            return pos0 - event.me_end - 1
        if region == "ex":
            return pos0 - event.me_start + 2
        if region == "c1":
            return pos0 - event.c1
        if region == "c2":
            return pos0 - event.c2
    else:
        if region == "up":
            return event.me_end - pos0
        if region == "dn":
            return event.me_start - pos0 - 1
        if region == "ex":
            return event.me_end - pos0 + 2
        if region == "c1":
            return pos0 - event.c1
        if region == "c2":
            return pos0 - event.c2
    raise ValueError(f"Unknown region: {region}")


def _query_hits_for_allele(
    data: PredictionData,
    record: pysam.VariantRecord,
    alt: str,
) -> list[QueryHit]:
    chrom = normalize_chrom(record.contig)
    pos = int(record.pos)
    ref = record.ref.upper()
    alt = alt.upper()
    variant = f"{chrom}:{pos}:{ref}:{alt}"
    variant_type = classify_variant(ref, alt)

    start0 = pos - 1
    end0 = start0 + max(1, len(ref)) - 1
    overlaps = {
        (interval.event, interval.region): interval
        for interval in data.find_overlaps(chrom, start0, end0)
    }

    exact_indels = data.indels_by_variant.get(variant, []) if variant_type == "Indel" else []
    for row in exact_indels:
        key = (row["event"], row["var_MIC_region"])
        if key not in overlaps and row["event"] in data.events:
            overlaps[key] = Interval(row["event"], row["var_MIC_region"], start0, end0)

    hits: list[QueryHit] = []
    indel_predictions = {
        (row["event"], row["var_MIC_region"]): row for row in exact_indels
    }
    for event_id, region in sorted(overlaps):
        event = data.events[event_id]
        hit = QueryHit(
            chrom=chrom,
            pos=pos,
            record_id=record.id or ".",
            ref=ref,
            alt=alt,
            variant=variant,
            event=event_id,
            region=region,
            relative_position=code_model_relative_position(event, region, start0),
            mic_distance=st12_mic_distance(event, region, start0),
            variant_type=variant_type,
        )
        indel_row = indel_predictions.get((event_id, region))
        if indel_row is not None:
            hit.prediction = {
                "mt_code_score": indel_row["mt_code_score"],
                "wt_code_score": indel_row["wt_code_score"],
                "logit_scaling_factor": indel_row[data.scaling_factor_column],
                "delta_code_model": indel_row[data.delta_code_model_column],
            }
            hit.mic_distance = indel_row["MIC_dist"]
            hit.prediction_status = "FOUND"
        hits.append(hit)
    return hits


def read_queries(
    input_vcf: Path,
    data: PredictionData,
    reference: ReferenceValidator | None,
    stats: Counter,
) -> tuple[list[QueryHit], set[tuple[str, str, int, str, str]]]:
    hits: list[QueryHit] = []
    wanted_snv: set[tuple[str, str, int, str, str]] = set()
    warned_reference: set[str] = set()

    with pysam.VariantFile(str(input_vcf)) as variant_file:
        for record in variant_file:
            stats["input_records"] += 1
            for alt_value in record.alts or ():
                stats["input_alt_alleles"] += 1
                alt = str(alt_value).upper()
                allele_hits = _query_hits_for_allele(data, record, alt)
                if not allele_hits:
                    stats["omitted_outside_coverage"] += 1
                    continue

                ref = record.ref.upper()
                unsupported = classify_variant(ref, alt) == "Unsupported"
                reference_state = "NOT_REQUESTED"
                reference_message = None
                if reference is not None and not unsupported:
                    reference_state, reference_message = reference.check(record.contig, record.pos, ref)

                variant = allele_hits[0].variant
                if unsupported:
                    for hit in allele_hits:
                        hit.prediction = None
                        hit.prediction_status = "ERROR_UNSUPPORTED_ALLELE"
                        hit.reference_check = "NOT_RUN"
                elif reference_state != "PASS" and reference is not None:
                    status = (
                        "ERROR_REF_MISMATCH"
                        if reference_state == "FAIL"
                        else "ERROR_REFERENCE_LOOKUP"
                    )
                    for hit in allele_hits:
                        hit.prediction = None
                        hit.prediction_status = status
                        hit.reference_check = reference_state
                    if variant not in warned_reference:
                        print(
                            f"WARNING: {variant}: {reference_message}",
                            file=sys.stderr,
                        )
                        warned_reference.add(variant)
                else:
                    for hit in allele_hits:
                        hit.reference_check = reference_state
                        if hit.variant_type == "SNV":
                            wanted_snv.add(
                                (hit.event, hit.region, hit.relative_position, hit.ref, hit.alt)
                            )
                hits.extend(allele_hits)

    return hits, wanted_snv


def attach_snv_predictions(
    hits: list[QueryHit],
    data: PredictionData,
    snv_predictions: dict[tuple[str, str, int, str, str], str],
) -> None:
    for hit in hits:
        if hit.variant_type != "SNV" or hit.prediction_status.startswith("ERROR"):
            continue
        key = (hit.event, hit.region, hit.relative_position, hit.ref, hit.alt)
        mt_score_text = snv_predictions.get(key)
        wt_score_text = data.wt_by_event.get(hit.event)
        scale = data.scale_by_event.get(hit.event)
        if mt_score_text is None or wt_score_text is None or scale in {None, 0}:
            continue
        mt_score = float(mt_score_text)
        wt_score = float(wt_score_text)
        delta = (logit2(mt_score) - logit2(wt_score)) / float(scale)
        hit.prediction = {
            "mt_code_score": mt_score_text,
            "wt_code_score": wt_score_text,
            "logit_scaling_factor": str(scale),
            "delta_code_model": str(delta),
        }
        hit.prediction_status = "FOUND"


def _open_output(path: Path) -> tuple[TextIO, bool]:
    if str(path) == "-":
        return sys.stdout, False
    path.parent.mkdir(parents=True, exist_ok=True)
    return path.open("w", newline=""), True


def write_results(output_path: Path, hits: Iterable[QueryHit], data: PredictionData, stats: Counter) -> None:
    handle, should_close = _open_output(output_path)
    try:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for hit in hits:
            event = data.events[hit.event]
            if hit.prediction_status.startswith("ERROR"):
                scores = {
                    "mt_code_score": "ERROR",
                    "wt_code_score": "ERROR",
                    "logit_scaling_factor": "ERROR",
                    "delta_code_model": "ERROR",
                }
                stats["errors"] += 1
            elif hit.prediction is None:
                scores = {
                    "mt_code_score": "NA",
                    "wt_code_score": "NA",
                    "logit_scaling_factor": (
                        str(data.scale_by_event[hit.event])
                        if hit.variant_type == "SNV"
                        else "NA"
                    ),
                    "delta_code_model": "NA",
                }
                stats["na_predictions"] += 1
            else:
                scores = hit.prediction
                stats["found_predictions"] += 1

            writer.writerow(
                {
                    "variant": hit.variant,
                    "geneName": event.gene_name,
                    "var_type": hit.variant_type,
                    "event": hit.event,
                    "ME-ID": event.microexon_id,
                    "MIC_group": event.microexon_group,
                    "MIC_dist": hit.mic_distance,
                    "var_MIC_region": hit.region,
                    **scores,
                    "assembly": "hg38",
                    "VCF_ID": hit.record_id,
                    "prediction_status": hit.prediction_status,
                    "ref_ok": {
                        "NOT_REQUESTED": "NA",
                        "PASS": "True",
                        "FAIL": "False",
                    }.get(hit.reference_check, hit.reference_check),
                }
            )
            stats["output_rows"] += 1
    finally:
        if should_close:
            handle.close()


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Retrieve precomputed microexon code-model predictions for variants in a "
            "GRCh38 VCF, compressed VCF, or BCF file."
        )
    )
    parser.add_argument(
        "-i",
        "--input-vcf",
        required=True,
        type=Path,
        help="Input GRCh38 VCF, VCF.GZ, or BCF file.",
    )
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        type=Path,
        help="Output TSV path. Use '-' to write to stdout.",
    )
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=Path(__file__).resolve().parent / "data",
        help="Prediction data directory (default: data/ beside this script).",
    )
    parser.add_argument(
        "--reference-fasta",
        type=Path,
        help="Optional indexed hg38 FASTA used to validate VCF REF alleles.",
    )
    parser.add_argument(
        "--delta-scaling",
        choices=tuple(DELTA_SCALING_COLUMNS),
        default="all-variants",
        help=(
            "Delta-score scaling to return: 'all-variants' (default) uses factors "
            "recalculated across the expanded variant set; 'main-analysis' uses the "
            "factors from the paper's main analyses."
        ),
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    stats: Counter = Counter()
    reference: ReferenceValidator | None = None
    try:
        data = PredictionData(args.data_dir, delta_scaling=args.delta_scaling)
        if args.reference_fasta is not None:
            reference = ReferenceValidator(args.reference_fasta)
        hits, wanted_snv = read_queries(args.input_vcf, data, reference, stats)
        snv_predictions = data.stream_snv_predictions(wanted_snv)
        attach_snv_predictions(hits, data, snv_predictions)
        write_results(args.output, hits, data, stats)
    except (FileNotFoundError, ValueError, OSError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    finally:
        if reference is not None:
            reference.close()

    print(
        "Summary: "
        f"records={stats['input_records']:,}; "
        f"ALT_alleles={stats['input_alt_alleles']:,}; "
        f"output_rows={stats['output_rows']:,}; "
        f"found={stats['found_predictions']:,}; "
        f"omitted_outside_coverage={stats['omitted_outside_coverage']:,}; "
        f"NA={stats['na_predictions']:,}; "
        f"errors={stats['errors']:,}",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
