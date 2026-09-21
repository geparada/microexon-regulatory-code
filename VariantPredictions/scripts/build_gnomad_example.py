#!/usr/bin/env python3
"""Build a diverse 100-variant example from public gnomAD v4.1 joint sites."""

from __future__ import annotations

import argparse
import sys
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

import pysam

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR.parent))

from variant_predictions import (  # noqa: E402
    PredictionData,
    classify_variant,
    code_model_relative_position,
    logit2,
    normalize_chrom,
)


GNOMAD_URL = (
    "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/joint/"
    "gnomad.joint.v4.1.sites.{chrom}.vcf.bgz"
)
REGION_ORDER = ("up", "dn", "c1", "c2", "ex")
FLANK_REGIONS = {"up", "dn"}


@dataclass
class CandidateMatch:
    event: str
    region: str
    relative_position: int
    microexon_group: str
    delta_code_model: float | None = None


@dataclass
class Candidate:
    chrom: str
    pos: int
    record_id: str
    ref: str
    alt: str
    variant_type: str
    af: str
    ac: str
    an: str
    matches: tuple[CandidateMatch, ...]
    source_url: str
    has_precomputed_indel: bool

    @property
    def variant(self) -> str:
        return f"{self.chrom}:{self.pos}:{self.ref}:{self.alt}"

    @property
    def key(self) -> tuple[int, int, str, str]:
        bare_chrom = self.chrom.removeprefix("chr")
        chrom_number = int(bare_chrom) if bare_chrom.isdigit() else 99
        return chrom_number, self.pos, self.ref, self.alt

    @property
    def events(self) -> set[str]:
        return {match.event for match in self.matches}

    @property
    def regions(self) -> set[str]:
        return {match.region for match in self.matches}

    @property
    def novel_events(self) -> set[str]:
        return {
            match.event
            for match in self.matches
            if match.microexon_group == "new_mic"
        }

    @property
    def best_flank_delta(self) -> float | None:
        values = [
            match.delta_code_model
            for match in self.matches
            if match.region in FLANK_REGIONS and match.delta_code_model is not None
        ]
        return min(values) if values else None


def natural_chromosomes() -> list[str]:
    return [f"chr{i}" for i in range(1, 23)] + ["chrX"]


def event_intervals(data: PredictionData, chrom: str) -> list:
    intervals = {
        (interval.event, interval.region): interval
        for values in data.interval_bins.get(chrom, {}).values()
        for interval in values
    }
    region_rank = {region: rank for rank, region in enumerate(REGION_ORDER)}
    return sorted(
        intervals.values(),
        key=lambda interval: (interval.event, region_rank[interval.region]),
    )


def alt_info_value(
    record: pysam.VariantRecord, names: tuple[str, ...], alt_index: int
) -> str:
    for name in names:
        if name not in record.info:
            continue
        value = record.info[name]
        if isinstance(value, tuple):
            if alt_index < len(value):
                return str(value[alt_index])
            continue
        return str(value)
    return "."


def is_pass(record: pysam.VariantRecord) -> bool:
    filters = set(record.filter.keys())
    return not filters or filters == {"PASS"}


def make_candidate(
    data: PredictionData,
    record: pysam.VariantRecord,
    alt_index: int,
    alt_value: str,
    source_url: str,
) -> Candidate | None:
    ref = record.ref.upper()
    alt = str(alt_value).upper()
    variant_type = classify_variant(ref, alt)
    if variant_type == "Unsupported":
        return None

    chrom = normalize_chrom(record.contig)
    start0 = record.pos - 1
    overlaps = data.find_overlaps(chrom, start0, start0 + len(ref) - 1)
    if not overlaps:
        return None

    matches = tuple(
        CandidateMatch(
            event=overlap.event,
            region=overlap.region,
            relative_position=code_model_relative_position(
                data.events[overlap.event], overlap.region, start0
            ),
            microexon_group=data.events[overlap.event].microexon_group,
        )
        for overlap in overlaps
    )
    variant = f"{chrom}:{record.pos}:{ref}:{alt}"
    return Candidate(
        chrom=chrom,
        pos=record.pos,
        record_id=record.id or ".",
        ref=ref,
        alt=alt,
        variant_type=variant_type,
        af=alt_info_value(record, ("AF_joint", "AF"), alt_index),
        ac=alt_info_value(record, ("AC_joint", "AC"), alt_index),
        an=alt_info_value(record, ("AN_joint", "AN"), alt_index),
        matches=matches,
        source_url=source_url,
        has_precomputed_indel=(
            variant_type == "Indel" and variant in data.indels_by_variant
        ),
    )


def match_distance(candidate: Candidate, event: str, region: str) -> int:
    return min(
        abs(match.relative_position)
        for match in candidate.matches
        if match.event == event and match.region == region
    )


def collect_candidates(
    data: PredictionData,
    target_pool: int,
    target_events: int,
) -> list[Candidate]:
    candidates: dict[tuple[str, int, str, str], Candidate] = {}
    queried_chromosomes = set()
    region_counts: Counter = Counter()

    for chrom in natural_chromosomes():
        intervals = event_intervals(data, chrom)
        if not intervals:
            continue
        source_url = GNOMAD_URL.format(chrom=chrom)
        print(f"Querying {chrom} from gnomAD v4.1 joint", file=sys.stderr)
        with pysam.VariantFile(source_url) as variant_file:
            source_contig = (
                chrom
                if chrom in variant_file.header.contigs
                else chrom.removeprefix("chr")
            )
            for interval in intervals:
                if (
                    interval.region not in FLANK_REGIONS
                    and region_counts[interval.region] >= 50
                ):
                    continue

                local_candidates: dict[
                    tuple[str, int, str, str], Candidate
                ] = {}
                try:
                    records = variant_file.fetch(
                        source_contig, interval.start0, interval.end0 + 1
                    )
                except (ValueError, OSError) as exc:
                    raise RuntimeError(f"Could not query {source_url}: {exc}") from exc

                for record in records:
                    if not is_pass(record):
                        continue
                    for alt_index, alt_value in enumerate(record.alts or ()):
                        candidate = make_candidate(
                            data, record, alt_index, alt_value, source_url
                        )
                        if candidate is None:
                            continue
                        if not any(
                            match.event == interval.event
                            and match.region == interval.region
                            for match in candidate.matches
                        ):
                            continue
                        key = (
                            candidate.chrom,
                            candidate.pos,
                            candidate.ref,
                            candidate.alt,
                        )
                        local_candidates[key] = candidate

                snvs = sorted(
                    (
                        candidate
                        for candidate in local_candidates.values()
                        if candidate.variant_type == "SNV"
                    ),
                    key=lambda candidate: (
                        match_distance(candidate, interval.event, interval.region),
                        candidate.key,
                    ),
                )[:12]
                indels = sorted(
                    (
                        candidate
                        for candidate in local_candidates.values()
                        if candidate.variant_type == "Indel"
                    ),
                    key=lambda candidate: (
                        not candidate.has_precomputed_indel,
                        match_distance(candidate, interval.event, interval.region),
                        candidate.key,
                    ),
                )[:3]

                for candidate in snvs + indels:
                    key = (
                        candidate.chrom,
                        candidate.pos,
                        candidate.ref,
                        candidate.alt,
                    )
                    if key in candidates:
                        continue
                    candidates[key] = candidate
                    region_counts.update(candidate.regions)

        queried_chromosomes.add(chrom)
        represented_events = {
            event for candidate in candidates.values() for event in candidate.events
        }
        novel_events = {
            event for candidate in candidates.values() for event in candidate.novel_events
        }
        indel_count = sum(
            candidate.variant_type == "Indel" for candidate in candidates.values()
        )
        if (
            len(candidates) >= target_pool
            and len(represented_events) >= target_events
            and len(novel_events) >= 50
            and indel_count >= 30
            and len(queried_chromosomes) >= 5
            and all(region_counts[region] >= 15 for region in REGION_ORDER)
        ):
            break

    return sorted(candidates.values(), key=lambda candidate: candidate.key)


def attach_snv_scores(data: PredictionData, candidates: list[Candidate]) -> None:
    wanted = {
        (
            match.event,
            match.region,
            match.relative_position,
            candidate.ref,
            candidate.alt,
        )
        for candidate in candidates
        if candidate.variant_type == "SNV"
        for match in candidate.matches
    }
    predictions = data.stream_snv_predictions(wanted)
    for candidate in candidates:
        if candidate.variant_type != "SNV":
            continue
        for match in candidate.matches:
            key = (
                match.event,
                match.region,
                match.relative_position,
                candidate.ref,
                candidate.alt,
            )
            mt_score = predictions.get(key)
            if mt_score is None:
                continue
            wt_score = float(data.wt_by_event[match.event])
            scale = data.scale_by_event[match.event]
            match.delta_code_model = (
                logit2(float(mt_score)) - logit2(wt_score)
            ) / scale


def score_sort_key(candidate: Candidate) -> tuple:
    delta = candidate.best_flank_delta
    return (delta is None, delta if delta is not None else 0, candidate.key)


def select_candidates(candidates: list[Candidate], count: int) -> list[Candidate]:
    if len(candidates) < count:
        raise ValueError(f"Only {len(candidates)} eligible gnomAD variants were found")

    selected: list[Candidate] = []
    selected_keys = set()
    covered_events = set()

    def add(candidate: Candidate, require_new_event: bool = True) -> bool:
        key = (candidate.chrom, candidate.pos, candidate.ref, candidate.alt)
        if key in selected_keys or len(selected) >= count:
            return False
        if require_new_event and not (candidate.events - covered_events):
            return False
        selected.append(candidate)
        selected_keys.add(key)
        covered_events.update(candidate.events)
        return True

    def add_until(pool: list[Candidate], metric, target: int) -> None:
        for require_new_event in (True, False):
            for candidate in pool:
                if metric() >= target or len(selected) >= count:
                    return
                add(candidate, require_new_event=require_new_event)

    flank_snvs = sorted(
        (
            candidate
            for candidate in candidates
            if candidate.variant_type == "SNV"
            and candidate.regions & FLANK_REGIONS
            and candidate.best_flank_delta is not None
        ),
        key=score_sort_key,
    )
    strongly_negative = [
        candidate
        for candidate in flank_snvs
        if candidate.best_flank_delta is not None
        and candidate.best_flank_delta <= -0.5
    ]
    add_until(
        strongly_negative,
        lambda: sum(
            candidate.best_flank_delta is not None
            and candidate.best_flank_delta <= -0.5
            for candidate in selected
        ),
        15,
    )

    indels = sorted(
        (candidate for candidate in candidates if candidate.variant_type == "Indel"),
        key=lambda candidate: (
            not candidate.has_precomputed_indel,
            not bool(candidate.regions & FLANK_REGIONS),
            not bool(candidate.novel_events),
            candidate.key,
        ),
    )
    add_until(
        indels,
        lambda: sum(candidate.variant_type == "Indel" for candidate in selected),
        10,
    )

    novel_flanks = sorted(
        (
            candidate for candidate in flank_snvs if candidate.novel_events
        ),
        key=score_sort_key,
    )
    add_until(
        novel_flanks,
        lambda: len(
            {
                event
                for candidate in selected
                for event in candidate.novel_events
            }
        ),
        20,
    )

    for region in ("c1", "c2", "ex"):
        region_pool = sorted(
            (candidate for candidate in candidates if region in candidate.regions),
            key=score_sort_key,
        )
        add_until(
            region_pool,
            lambda region=region: sum(
                region in candidate.regions for candidate in selected
            ),
            3,
        )

    add_until(
        flank_snvs,
        lambda: sum(
            bool(candidate.regions & FLANK_REGIONS) for candidate in selected
        ),
        80,
    )

    remaining = sorted(candidates, key=score_sort_key)
    add_until(remaining, lambda: len(selected), count)

    selected = sorted(selected, key=lambda candidate: candidate.key)
    selected_events = {
        event for candidate in selected for event in candidate.events
    }
    selected_novel_events = {
        event for candidate in selected for event in candidate.novel_events
    }
    selected_regions = {
        region for candidate in selected for region in candidate.regions
    }
    strong_negative_count = sum(
        candidate.best_flank_delta is not None
        and candidate.best_flank_delta <= -0.5
        for candidate in selected
    )
    flank_count = sum(
        bool(candidate.regions & FLANK_REGIONS) for candidate in selected
    )
    indel_count = sum(candidate.variant_type == "Indel" for candidate in selected)

    if len(selected) != count:
        raise ValueError(f"Could not select exactly {count} variants")
    if len(selected_events) < 80:
        raise ValueError(
            f"Only {len(selected_events)} microexon events were represented"
        )
    if len(selected_novel_events) < 15:
        raise ValueError(
            f"Only {len(selected_novel_events)} novel microexons were represented"
        )
    if flank_count < 75:
        raise ValueError(f"Only {flank_count} selected variants overlap up/dn regions")
    if indel_count < 10:
        raise ValueError(f"Only {indel_count} indels were selected")
    if strong_negative_count < 5:
        raise ValueError(
            f"Only {strong_negative_count} variants have delta_code_model <= -0.5"
        )
    if selected_regions != set(REGION_ORDER):
        raise ValueError(f"Selected region labels are incomplete: {selected_regions}")
    return selected


def chrom_sort_key(chrom: str) -> int:
    bare_chrom = chrom.removeprefix("chr")
    return int(bare_chrom) if bare_chrom.isdigit() else 99


def write_vcf(path: Path, candidates: list[Candidate]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    source_urls = sorted({candidate.source_url for candidate in candidates})
    contigs = sorted(
        {candidate.chrom for candidate in candidates}, key=chrom_sort_key
    )
    with path.open("w", newline="") as handle:
        handle.write("##fileformat=VCFv4.2\n")
        handle.write("##reference=GRCh38\n")
        handle.write("##source=gnomAD_v4.1_joint_microexon_example\n")
        for url in source_urls:
            handle.write(f"##gnomADSource={url}\n")
        for contig in contigs:
            handle.write(f"##contig=<ID={contig}>\n")
        handle.write(
            '##INFO=<ID=AF,Number=A,Type=Float,Description="gnomAD v4.1 joint allele frequency">\n'
        )
        handle.write(
            '##INFO=<ID=AC,Number=A,Type=Integer,Description="gnomAD v4.1 joint allele count">\n'
        )
        handle.write(
            '##INFO=<ID=AN,Number=1,Type=Integer,Description="gnomAD v4.1 joint allele number">\n'
        )
        handle.write(
            '##INFO=<ID=MIC_EVENT,Number=.,Type=String,Description="Overlapping microexon event identifiers">\n'
        )
        handle.write(
            '##INFO=<ID=MIC_GROUP,Number=.,Type=String,Description="Overlapping microexon groups">\n'
        )
        handle.write(
            '##INFO=<ID=MIC_REGION,Number=.,Type=String,Description="Overlapping microexon region labels">\n'
        )
        handle.write(
            '##INFO=<ID=MIN_DELTA_CODE_MODEL,Number=1,Type=Float,Description="Most negative available up/dn delta code-model score">\n'
        )
        handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for candidate in candidates:
            events = sorted(candidate.events)
            groups = sorted(
                {match.microexon_group for match in candidate.matches}
            )
            regions = [
                region for region in REGION_ORDER if region in candidate.regions
            ]
            delta = (
                str(candidate.best_flank_delta)
                if candidate.best_flank_delta is not None
                else "."
            )
            info = (
                f"AF={candidate.af};AC={candidate.ac};AN={candidate.an};"
                f"MIC_EVENT={','.join(events)};MIC_GROUP={','.join(groups)};"
                f"MIC_REGION={','.join(regions)};MIN_DELTA_CODE_MODEL={delta}"
            )
            handle.write(
                f"{candidate.chrom}\t{candidate.pos}\t{candidate.record_id}\t"
                f"{candidate.ref}\t{candidate.alt}\t.\tPASS\t{info}\n"
            )


def describe_selection(candidates: list[Candidate]) -> str:
    events = {event for candidate in candidates for event in candidate.events}
    novel_events = {
        event for candidate in candidates for event in candidate.novel_events
    }
    flank_count = sum(
        bool(candidate.regions & FLANK_REGIONS) for candidate in candidates
    )
    negative_count = sum(
        candidate.best_flank_delta is not None
        and candidate.best_flank_delta <= -0.5
        for candidate in candidates
    )
    return (
        f"variants={len(candidates)}, events={len(events)}, "
        f"novel_events={len(novel_events)}, up_or_dn={flank_count}, "
        f"SNVs={sum(x.variant_type == 'SNV' for x in candidates)}, "
        f"indels={sum(x.variant_type == 'Indel' for x in candidates)}, "
        f"delta_le_-0.5={negative_count}"
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--count", type=int, default=100)
    parser.add_argument("--candidate-pool", type=int, default=1_200)
    parser.add_argument("--target-events", type=int, default=160)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    data = PredictionData(args.data_dir)
    candidates = collect_candidates(
        data,
        target_pool=max(args.candidate_pool, args.count),
        target_events=max(args.target_events, 80),
    )
    attach_snv_scores(data, candidates)
    print(f"Candidate pool: {describe_selection(candidates)}", file=sys.stderr)
    selected = select_candidates(candidates, args.count)
    write_vcf(args.output, selected)
    print(f"Wrote {args.output}: {describe_selection(selected)}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
