from __future__ import annotations

import argparse
import csv
import gzip
import sys
import tempfile
import unittest
from pathlib import Path


TOOL_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(TOOL_DIR / "scripts"))

import build_prediction_data as build  # noqa: E402


class BuildPredictionDataTests(unittest.TestCase):
    def test_bundle_membership_and_complete_event_scaling(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            ism = root / "MIC_ISM.gz"
            wt = root / "MIC_ISM_wt.gz"
            st12 = root / "ST12.csv"
            metadata = root / "microexons.csv"
            output = root / "output"

            self._write_gzip_csv(
                ism,
                [
                    "event",
                    "region",
                    "rel_pos",
                    "ref",
                    "alt",
                    "mt_score",
                    "delta_logit_score",
                    "transformed_delta_scores",
                ],
                [
                    self._ism_row("ISM_A", "-2"),
                    self._ism_row("ISM_A", "1"),
                    self._ism_row("ISM_B", "-3"),
                ],
            )
            self._write_gzip_csv(
                wt,
                ["event", "wt_score"],
                [
                    {"event": "ISM_A", "wt_score": "0.8"},
                    {"event": "ISM_B", "wt_score": "0.7"},
                    {"event": "WT_ONLY", "wt_score": "0.6"},
                ],
            )
            self._write_csv(
                metadata,
                build.EVENT_COLUMNS,
                [
                    self._metadata_row("ISM_A"),
                    self._metadata_row("ISM_B"),
                    self._metadata_row("WT_ONLY"),
                    self._metadata_row("ST12_ONLY"),
                ],
            )
            self._write_csv(
                st12,
                [
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
                    "var_type",
                    "assembly",
                ],
                [
                    self._st12_row("ISM_A", "chr1:1001:A:AT", "2", "4"),
                    self._st12_row("ST12_ONLY", "chr1:1002:A:AT", "5", "6"),
                ],
            )

            manifest = build.build_data(
                argparse.Namespace(
                    mic_ism=ism,
                    mic_wt=wt,
                    st12=st12,
                    microexons=metadata,
                    output_dir=output,
                )
            )

            self.assertEqual(manifest["covered_events"], 2)
            self.assertEqual(manifest["files"]["MIC_indels.csv.gz"]["rows"], 1)
            self.assertEqual(
                self._read_events(output / "MIC_ISM_wt.gz"), {"ISM_A", "ISM_B"}
            )
            self.assertEqual(
                self._read_events(output / "microexon_events.hg38.csv.gz"),
                {"ISM_A", "ISM_B"},
            )

            with gzip.open(
                output / "MIC_event_scaling.csv.gz", "rt", newline=""
            ) as handle:
                scaling = {row["event"]: row for row in csv.DictReader(handle)}
            self.assertEqual(set(scaling), {"ISM_A", "ISM_B"})
            self.assertEqual(scaling["ISM_A"]["logit_scaling_factor"], "2.0")
            self.assertEqual(
                scaling["ISM_A"]["logit_scaling_factor_all_vars"], "4.0"
            )
            self.assertEqual(scaling["ISM_B"]["logit_scaling_factor"], "3.0")
            self.assertEqual(
                scaling["ISM_B"]["logit_scaling_factor_all_vars"], "3.0"
            )

    @staticmethod
    def _ism_row(event: str, delta: str) -> dict[str, str]:
        return {
            "event": event,
            "region": "c1",
            "rel_pos": "0",
            "ref": "A",
            "alt": "G",
            "mt_score": "0.5",
            "delta_logit_score": delta,
            "transformed_delta_scores": "0",
        }

    @staticmethod
    def _metadata_row(event: str) -> dict[str, str]:
        return {
            "event": event,
            "chrom": "chr1",
            "strand": "+",
            "upIntStart": "1000",
            "upIntEnd": "2000",
            "dnIntStart": "2003",
            "dnIntEnd": "3000",
            "geneName": event,
            "lengthDiff": "3",
            "group": "known_mic",
        }

    @staticmethod
    def _st12_row(
        event: str, variant: str, historical_factor: str, all_vars_factor: str
    ) -> dict[str, str]:
        return {
            "variant": variant,
            "event": event,
            "var_MIC_region": "c1",
            "MIC_dist": "0",
            "mt_code_score": "0.4",
            "wt_code_score": "0.8",
            "logit_scaling_factor": historical_factor,
            "delta_code_model": "-1.0",
            "logit_scaling_factor_all_vars": all_vars_factor,
            "delta_code_model_all_vars": "-0.646240625180289",
            "var_type": "Indel",
            "assembly": "hg38",
        }

    @staticmethod
    def _write_csv(path: Path, columns: list[str], rows: list[dict[str, str]]) -> None:
        with path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=columns, lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)

    @staticmethod
    def _write_gzip_csv(
        path: Path, columns: list[str], rows: list[dict[str, str]]
    ) -> None:
        with gzip.open(path, "wt", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=columns, lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)

    @staticmethod
    def _read_events(path: Path) -> set[str]:
        with gzip.open(path, "rt", newline="") as handle:
            return {row["event"] for row in csv.DictReader(handle)}


if __name__ == "__main__":
    unittest.main()
