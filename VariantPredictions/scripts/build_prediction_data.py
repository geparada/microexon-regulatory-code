#!/usr/bin/env python3
"""Build the compact runtime data bundle for variant_predictions.py."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import math
import shutil
from contextlib import contextmanager
from pathlib import Path
from typing import Iterator, TextIO


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


@contextmanager
def deterministic_gzip_text(path: Path) -> Iterator[TextIO]:
    with path.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as compressed:
            with io.TextIOWrapper(compressed, encoding="utf-8", newline="") as text:
                yield text


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def require_columns(actual: list[str] | None, expected: list[str], source: Path) -> None:
    missing = [column for column in expected if column not in (actual or [])]
    if missing:
        raise ValueError(f"{source} is missing required columns: {', '.join(missing)}")


def validate_scale(
    scales: dict[str, float], event: str, value: str, column: str
) -> None:
    if value in {"", "NA"}:
        raise ValueError(f"Missing {column} for {event}")
    scale = float(value)
    if not math.isfinite(scale) or scale <= 0:
        raise ValueError(f"Invalid {column} for {event}: {value}")
    if event in scales and not math.isclose(
        scales[event], scale, rel_tol=0, abs_tol=1e-12
    ):
        raise ValueError(f"Conflicting {column} values for {event}")
    scales[event] = scale


def logit2(probability: str) -> float:
    value = float(probability)
    epsilon = 1e-12
    value = min(max(value, epsilon), 1 - epsilon)
    return math.log2(value / (1 - value))


def validate_indel_scores(row: dict[str, str]) -> None:
    raw_delta = logit2(row["mt_code_score"]) - logit2(row["wt_code_score"])
    historical = raw_delta / float(row["logit_scaling_factor"])
    expected = {
        "delta_code_model": max(-1.0, min(1.0, historical)),
        "delta_code_model_all_vars": raw_delta
        / float(row["logit_scaling_factor_all_vars"]),
    }
    for column, expected_value in expected.items():
        observed = float(row[column])
        if not math.isclose(observed, expected_value, rel_tol=1e-11, abs_tol=1e-12):
            raise ValueError(
                f"Indel row {row['variant']} has an inconsistent {column}: "
                f"observed={observed}, expected={expected_value}"
            )


def build_data(args: argparse.Namespace) -> dict:
    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    covered_events = set()
    ism_scale_by_event: dict[str, float] = {}
    ism_rows = 0
    with gzip.open(args.mic_ism, "rt", newline="") as handle:
        reader = csv.DictReader(handle)
        require_columns(reader.fieldnames, ISM_COLUMNS, args.mic_ism)
        for row in reader:
            event = row["event"]
            delta_logit = abs(float(row["delta_logit_score"]))
            if not math.isfinite(delta_logit):
                raise ValueError(f"Non-finite delta_logit_score for {event}")
            covered_events.add(event)
            ism_scale_by_event[event] = max(
                delta_logit, ism_scale_by_event.get(event, 0.0)
            )
            ism_rows += 1
    if not covered_events:
        raise ValueError(f"No events found in {args.mic_ism}")

    mic_ism_output = output_dir / "MIC_ISM.gz"
    shutil.copy2(args.mic_ism, mic_ism_output)

    wt_by_event = {}
    with gzip.open(args.mic_wt, "rt", newline="") as handle:
        reader = csv.DictReader(handle)
        require_columns(reader.fieldnames, ["event", "wt_score"], args.mic_wt)
        for row in reader:
            event = row["event"]
            if event not in covered_events:
                continue
            if event in wt_by_event and wt_by_event[event] != row["wt_score"]:
                raise ValueError(f"Conflicting wild-type scores for {event}")
            wt_by_event[event] = row["wt_score"]
    if set(wt_by_event) != covered_events:
        missing = sorted(covered_events - set(wt_by_event))
        raise ValueError(
            "MIC_ISM_wt.gz does not contain every MIC_ISM event: "
            f"missing={len(missing)}"
        )

    mic_wt_output = output_dir / "MIC_ISM_wt.gz"
    with deterministic_gzip_text(mic_wt_output) as output_handle:
        writer = csv.DictWriter(
            output_handle, fieldnames=["event", "wt_score"], lineterminator="\n"
        )
        writer.writeheader()
        for event in sorted(covered_events):
            writer.writerow({"event": event, "wt_score": wt_by_event[event]})

    metadata_output = output_dir / "microexon_events.hg38.csv.gz"
    metadata_rows = 0
    metadata_events = set()
    with args.microexons.open(newline="") as input_handle:
        reader = csv.DictReader(input_handle)
        require_columns(reader.fieldnames, EVENT_COLUMNS, args.microexons)
        with deterministic_gzip_text(metadata_output) as output_handle:
            writer = csv.DictWriter(output_handle, fieldnames=EVENT_COLUMNS, lineterminator="\n")
            writer.writeheader()
            for row in reader:
                if row["event"] not in covered_events:
                    continue
                writer.writerow({column: row[column] for column in EVENT_COLUMNS})
                metadata_rows += 1
                metadata_events.add(row["event"])
    if metadata_events != covered_events or metadata_rows != len(covered_events):
        raise ValueError(
            "Filtered metadata does not contain exactly one row per MIC_ISM event: "
            f"rows={metadata_rows}, events={len(metadata_events)}, expected={len(covered_events)}"
        )

    indel_output = output_dir / "MIC_indels.csv.gz"
    indel_rows = 0
    indel_variants = set()
    indel_events = set()
    historical_scale_by_event: dict[str, float] = {}
    all_vars_scale_by_event: dict[str, float] = {}
    with args.st12.open(newline="") as input_handle:
        reader = csv.DictReader(input_handle)
        required = set(INDEL_COLUMNS) | {"var_type", "assembly"}
        missing = sorted(required - set(reader.fieldnames or []))
        if missing:
            raise ValueError(f"Supplementary Table 12 is missing columns: {missing}")
        with deterministic_gzip_text(indel_output) as output_handle:
            writer = csv.DictWriter(output_handle, fieldnames=INDEL_COLUMNS, lineterminator="\n")
            writer.writeheader()
            for row in reader:
                event = row["event"]
                if event not in covered_events:
                    continue
                if row["assembly"] != "hg38":
                    raise ValueError(f"Unexpected assembly for {row['variant']}: {row['assembly']}")
                scale_values = (
                    row["logit_scaling_factor"],
                    row["logit_scaling_factor_all_vars"],
                )
                if any(value in {"", "NA"} for value in scale_values):
                    if row["var_type"] == "Indel":
                        raise ValueError(
                            f"Indel row {row['variant']} has missing scaling values"
                        )
                    continue
                validate_scale(
                    historical_scale_by_event,
                    event,
                    row["logit_scaling_factor"],
                    "logit_scaling_factor",
                )
                validate_scale(
                    all_vars_scale_by_event,
                    event,
                    row["logit_scaling_factor_all_vars"],
                    "logit_scaling_factor_all_vars",
                )
                if row["var_type"] != "Indel":
                    continue
                missing_values = [column for column in INDEL_COLUMNS if row[column] in {"", "NA"}]
                if missing_values:
                    raise ValueError(
                        f"Indel row {row['variant']} has missing values: {missing_values}"
                    )
                validate_indel_scores(row)
                writer.writerow({column: row[column] for column in INDEL_COLUMNS})
                indel_rows += 1
                indel_variants.add(row["variant"])
                indel_events.add(event)

    if set(historical_scale_by_event) != set(all_vars_scale_by_event):
        raise ValueError(
            "Historical and all-variant Supplementary Table 12 scaling coverage differs"
        )

    scaling_output = output_dir / "MIC_event_scaling.csv.gz"
    with deterministic_gzip_text(scaling_output) as output_handle:
        writer = csv.DictWriter(
            output_handle, fieldnames=SCALING_COLUMNS, lineterminator="\n"
        )
        writer.writeheader()
        for event in sorted(covered_events):
            if event in historical_scale_by_event:
                historical = historical_scale_by_event[event]
                all_vars = all_vars_scale_by_event[event]
            else:
                historical = ism_scale_by_event[event]
                all_vars = historical
                if historical <= 0:
                    raise ValueError(
                        f"Cannot derive a positive ISM scaling factor for {event}"
                    )
            writer.writerow(
                {
                    "event": event,
                    "logit_scaling_factor": str(historical),
                    "logit_scaling_factor_all_vars": str(all_vars),
                }
            )

    descriptions = {
        "MIC_ISM.gz": {"rows": ism_rows, "events": len(covered_events)},
        "MIC_ISM_wt.gz": {"rows": len(wt_by_event), "events": len(covered_events)},
        "microexon_events.hg38.csv.gz": {
            "rows": metadata_rows,
            "events": len(metadata_events),
        },
        "MIC_indels.csv.gz": {
            "rows": indel_rows,
            "events": len(indel_events),
            "unique_variants": len(indel_variants),
        },
        "MIC_event_scaling.csv.gz": {
            "rows": len(covered_events),
            "events": len(covered_events),
        },
    }
    files = {}
    for name, description in descriptions.items():
        path = output_dir / name
        files[name] = {
            **description,
            "bytes": path.stat().st_size,
            "sha256": sha256(path),
        }

    manifest = {
        "format_version": 2,
        "assembly": "GRCh38/hg38",
        "covered_events": len(covered_events),
        "sources": {
            "mic_ism": args.mic_ism.name,
            "mic_wt": args.mic_wt.name,
            "supplementary_table_12": args.st12.name,
            "microexons": args.microexons.name,
        },
        "files": files,
    }
    manifest_path = output_dir / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    return manifest


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mic-ism", required=True, type=Path)
    parser.add_argument("--mic-wt", required=True, type=Path)
    parser.add_argument(
        "--supplementary-table-12",
        "--st12",
        dest="st12",
        required=True,
        type=Path,
        help="CSV file for Supplementary Table 12 of the paper.",
    )
    parser.add_argument("--microexons", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    manifest = build_data(args)
    print(json.dumps(manifest, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
