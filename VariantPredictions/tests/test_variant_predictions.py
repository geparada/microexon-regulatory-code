from __future__ import annotations

import csv
import gzip
import io
import math
import sys
import tempfile
import unittest
from contextlib import redirect_stderr
from pathlib import Path

import pysam


TOOL_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(TOOL_DIR))

import variant_predictions as vp  # noqa: E402


class VariantPredictionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.temporary_directory = tempfile.TemporaryDirectory()
        cls.root = Path(cls.temporary_directory.name)
        cls.data_dir = cls.root / "data"
        cls.data_dir.mkdir()
        cls._write_synthetic_data()
        cls.data = vp.PredictionData(cls.data_dir)

    @classmethod
    def tearDownClass(cls) -> None:
        cls.temporary_directory.cleanup()

    @classmethod
    def _write_gzip_csv(
        cls, name: str, fieldnames: list[str], rows: list[dict[str, str]]
    ) -> None:
        with gzip.open(cls.data_dir / name, "wt", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)

    @classmethod
    def _write_synthetic_data(cls) -> None:
        cls._write_gzip_csv(
            "MIC_ISM_wt.gz",
            ["event", "wt_score"],
            [
                {"event": "PLUS", "wt_score": "0.8"},
                {"event": "MINUS", "wt_score": "0.7"},
            ],
        )
        cls._write_gzip_csv(
            "MIC_ISM.gz",
            vp.ISM_COLUMNS,
            [
                {
                    "event": event,
                    "region": region,
                    "rel_pos": str(relative_position),
                    "ref": "A",
                    "alt": "G",
                    "mt_score": mutant_score,
                    "delta_logit_score": "0",
                    "transformed_delta_scores": "0",
                }
                for event, region, relative_position, mutant_score in (
                    ("PLUS", "c1", 0, "0.4"),
                    ("PLUS", "ex", 1, "0.5"),
                    ("MINUS", "c1", 0, "0.6"),
                    ("MINUS", "ex", 1, "0.3"),
                )
            ],
        )
        cls._write_gzip_csv(
            "MIC_event_scaling.csv.gz",
            [
                "event",
                "logit_scaling_factor",
                "logit_scaling_factor_all_vars",
            ],
            [
                {
                    "event": "PLUS",
                    "logit_scaling_factor": "2",
                    "logit_scaling_factor_all_vars": "4",
                },
                {
                    "event": "MINUS",
                    "logit_scaling_factor": "2.5",
                    "logit_scaling_factor_all_vars": "5",
                },
            ],
        )
        cls._write_gzip_csv(
            "microexon_events.hg38.csv.gz",
            vp.EVENT_COLUMNS,
            [
                {
                    "event": "PLUS",
                    "chrom": "chr1",
                    "strand": "+",
                    "upIntStart": "1000",
                    "upIntEnd": "2000",
                    "dnIntStart": "2003",
                    "dnIntEnd": "3000",
                    "geneName": "GENE_PLUS",
                    "lengthDiff": "3",
                    "group": "known_mic",
                },
                {
                    "event": "MINUS",
                    "chrom": "1",
                    "strand": "-",
                    "upIntStart": "6003",
                    "upIntEnd": "1000",
                    "dnIntStart": "7000",
                    "dnIntEnd": "6000",
                    "geneName": "GENE_MINUS",
                    "lengthDiff": "3",
                    "group": "novel_mic",
                },
            ],
        )
        cls._write_gzip_csv(
            "MIC_indels.csv.gz",
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
            ],
            [
                {
                    "variant": "chr1:1001:A:AT",
                    "event": "PLUS",
                    "var_MIC_region": "c1",
                    "MIC_dist": "0",
                    "mt_code_score": "0.6",
                    "wt_code_score": "0.8",
                    "logit_scaling_factor": "2",
                    "delta_code_model": "-0.5283208335737186",
                    "logit_scaling_factor_all_vars": "4",
                    "delta_code_model_all_vars": "-0.2641604167868593",
                },
                {
                    "variant": "chr1:1001:A:AT",
                    "event": "MINUS",
                    "var_MIC_region": "c1",
                    "MIC_dist": "0",
                    "mt_code_score": "0.5",
                    "wt_code_score": "0.7",
                    "logit_scaling_factor": "2.5",
                    "delta_code_model": "-0.4889082860886104",
                    "logit_scaling_factor_all_vars": "5",
                    "delta_code_model_all_vars": "-0.2444541430443052",
                },
            ],
        )

    def write_vcf(self, name: str, body: str, samples: bool = False) -> Path:
        path = self.root / name
        sample_headers = "\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" if samples else ""
        sample_columns = "\tFORMAT\tSAMPLE" if samples else ""
        path.write_text(
            "##fileformat=VCFv4.2\n"
            "##contig=<ID=chr1,length=10000>\n"
            "##contig=<ID=1,length=10000>"
            f"{sample_headers}\n"
            f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO{sample_columns}\n"
            f"{body}"
        )
        return path

    def run_tool(
        self,
        input_path: Path,
        output_name: str,
        reference_fasta: Path | None = None,
        delta_scaling: str | None = None,
    ) -> tuple[list[dict[str, str]], str, bytes]:
        output_path = self.root / output_name
        arguments = [
            "--input-vcf",
            str(input_path),
            "--output",
            str(output_path),
            "--data-dir",
            str(self.data_dir),
        ]
        if reference_fasta is not None:
            arguments.extend(["--reference-fasta", str(reference_fasta)])
        if delta_scaling is not None:
            arguments.extend(["--delta-scaling", delta_scaling])
        standard_error = io.StringIO()
        with redirect_stderr(standard_error):
            exit_code = vp.main(arguments)
        self.assertEqual(exit_code, 0, standard_error.getvalue())
        raw = output_path.read_bytes()
        with output_path.open(newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        return rows, standard_error.getvalue(), raw

    def test_default_uses_all_variants_scaling(self) -> None:
        vcf = self.write_vcf(
            "default_scaling.vcf",
            "chr1\t1001\tsnv\tA\tG\t.\tPASS\t.\n"
            "chr1\t1001\tindel\tA\tAT\t.\tPASS\t.\n",
        )
        rows, _, _ = self.run_tool(vcf, "default_scaling.tsv")
        plus_snv = next(
            row for row in rows if row["VCF_ID"] == "snv" and row["event"] == "PLUS"
        )
        plus_indel = next(
            row for row in rows if row["VCF_ID"] == "indel" and row["event"] == "PLUS"
        )
        self.assertEqual(plus_snv["logit_scaling_factor"], "4.0")
        self.assertTrue(
            math.isclose(float(plus_snv["delta_code_model"]), -0.6462406251802891)
        )
        self.assertEqual(plus_indel["logit_scaling_factor"], "4")
        self.assertEqual(plus_indel["delta_code_model"], "-0.2641604167868593")

    def test_main_analysis_scaling_selects_historical_factor_and_score(self) -> None:
        vcf = self.write_vcf(
            "main_analysis_scaling.vcf",
            "chr1\t1001\tsnv\tA\tG\t.\tPASS\t.\n"
            "chr1\t1001\tindel\tA\tAT\t.\tPASS\t.\n",
        )
        rows, _, _ = self.run_tool(
            vcf, "main_analysis_scaling.tsv", delta_scaling="main-analysis"
        )
        plus_snv = next(
            row for row in rows if row["VCF_ID"] == "snv" and row["event"] == "PLUS"
        )
        plus_indel = next(
            row for row in rows if row["VCF_ID"] == "indel" and row["event"] == "PLUS"
        )
        self.assertEqual(plus_snv["logit_scaling_factor"], "2.0")
        self.assertTrue(
            math.isclose(float(plus_snv["delta_code_model"]), -1.2924812503605782)
        )
        self.assertGreater(abs(float(plus_snv["delta_code_model"])), 1)
        self.assertEqual(plus_indel["logit_scaling_factor"], "2")
        self.assertEqual(plus_indel["delta_code_model"], "-0.5283208335737186")

    def test_plain_bgzf_bcf_samples_and_multiallelic_records(self) -> None:
        plain = self.write_vcf(
            "formats.vcf",
            "1\t1001\tmulti\tA\tG,C\t.\tPASS\t.\tGT\t1/2\n",
            samples=True,
        )
        bgzf = self.root / "formats.vcf.gz"
        pysam.tabix_compress(str(plain), str(bgzf), force=True)
        pysam.tabix_index(str(bgzf), preset="vcf", force=True)

        bcf = self.root / "formats.bcf"
        with pysam.VariantFile(str(plain)) as source:
            with pysam.VariantFile(str(bcf), "wb", header=source.header) as destination:
                for record in source:
                    destination.write(record)

        baseline = None
        for input_path in (plain, bgzf, bcf):
            with self.subTest(suffix=input_path.suffix):
                rows, summary, raw = self.run_tool(
                    input_path, f"{input_path.name}.predictions.tsv"
                )
                self.assertEqual(len(rows), 4)
                self.assertEqual(
                    [row["prediction_status"] for row in rows],
                    ["FOUND", "FOUND", "NA_NOT_PRECOMPUTED", "NA_NOT_PRECOMPUTED"],
                )
                self.assertTrue(all(row["variant"].startswith("chr1:") for row in rows))
                self.assertIn("records=1", summary)
                self.assertIn("ALT_alleles=2", summary)
                if baseline is None:
                    baseline = raw
                else:
                    self.assertEqual(raw, baseline)

    def test_strands_region_boundaries_and_relative_positions(self) -> None:
        expected = {
            "PLUS": {
                "c1": (700, 1299),
                "up": (1700, 1999),
                "ex": (2000, 2002),
                "dn": (2003, 2302),
                "c2": (2700, 3299),
            },
            "MINUS": {
                "c1": (700, 1299),
                "dn": (5700, 5999),
                "ex": (6000, 6002),
                "up": (6003, 6302),
                "c2": (6700, 7299),
            },
        }
        for event_id, regions in expected.items():
            for region, (start0, end0) in regions.items():
                for pos0 in (start0, end0):
                    overlaps = {
                        (hit.event, hit.region)
                        for hit in self.data.find_overlaps("1", pos0, pos0)
                    }
                    self.assertIn((event_id, region), overlaps)
                outside = {
                    (hit.event, hit.region)
                    for hit in self.data.find_overlaps("chr1", end0 + 1, end0 + 1)
                }
                self.assertNotIn((event_id, region), outside)

        self.assertEqual(
            vp.code_model_relative_position(self.data.events["PLUS"], "ex", 2000), 1
        )
        self.assertEqual(
            vp.code_model_relative_position(self.data.events["MINUS"], "ex", 6002), 1
        )
        self.assertEqual(
            vp.code_model_relative_position(self.data.events["PLUS"], "c1", 1000), 0
        )
        self.assertEqual(
            vp.code_model_relative_position(self.data.events["MINUS"], "c1", 1000), 0
        )
        self.assertEqual(
            vp.st12_mic_distance(self.data.events["PLUS"], "dn", 2003), 0
        )
        self.assertEqual(
            vp.st12_mic_distance(self.data.events["PLUS"], "ex", 2000), 2
        )
        self.assertEqual(
            vp.st12_mic_distance(self.data.events["MINUS"], "c1", 999), -1
        )

    def test_one_to_many_indel_exact_match_na_unsupported_and_omission(self) -> None:
        vcf = self.write_vcf(
            "behaviors.vcf",
            "chr1\t1001\tindel\tA\tAT\t.\tPASS\t.\n"
            "chr1\t1002\tmissing_indel\tA\tAT\t.\tPASS\t.\n"
            "chr1\t1001\tmissing\tA\tC\t.\tPASS\t.\n"
            "chr1\t1001\tsymbolic\tA\t<DEL>\t.\tPASS\t.\n"
            "chr1\t9001\toutside\tA\tG\t.\tPASS\t.\n",
        )
        rows, summary, _ = self.run_tool(vcf, "behaviors.tsv")
        by_id: dict[str, list[dict[str, str]]] = {}
        for row in rows:
            by_id.setdefault(row["VCF_ID"], []).append(row)

        self.assertEqual(len(by_id["indel"]), 2)
        self.assertEqual({row["event"] for row in by_id["indel"]}, {"PLUS", "MINUS"})
        self.assertTrue(all(row["prediction_status"] == "FOUND" for row in by_id["indel"]))
        self.assertEqual(len(by_id["missing_indel"]), 2)
        self.assertTrue(
            all(
                row["prediction_status"] == "NA_NOT_PRECOMPUTED"
                for row in by_id["missing_indel"]
            )
        )
        self.assertEqual(len(by_id["missing"]), 2)
        self.assertTrue(
            all(row["prediction_status"] == "NA_NOT_PRECOMPUTED" for row in by_id["missing"])
        )
        self.assertEqual(len(by_id["symbolic"]), 2)
        self.assertTrue(
            all(row["prediction_status"] == "ERROR_UNSUPPORTED_ALLELE" for row in by_id["symbolic"])
        )
        self.assertTrue(all(row["wt_code_score"] == "ERROR" for row in by_id["symbolic"]))
        self.assertNotIn("outside", by_id)
        self.assertIn("omitted_outside_coverage=1", summary)

    def test_reference_validation_pass_and_mismatch(self) -> None:
        fasta = self.root / "reference.fa"
        fasta.write_text(">chr1\n" + "A" * 10000 + "\n")
        pysam.faidx(str(fasta))
        vcf = self.write_vcf(
            "reference.vcf",
            "chr1\t1001\tpass\tA\tG\t.\tPASS\t.\n"
            "chr1\t1001\tmismatch\tC\tG\t.\tPASS\t.\n",
        )
        rows, standard_error, _ = self.run_tool(vcf, "reference.tsv", fasta)
        passing = [row for row in rows if row["VCF_ID"] == "pass"]
        mismatch = [row for row in rows if row["VCF_ID"] == "mismatch"]
        self.assertTrue(all(row["ref_ok"] == "True" for row in passing))
        self.assertTrue(all(row["prediction_status"] == "FOUND" for row in passing))
        self.assertTrue(all(row["ref_ok"] == "False" for row in mismatch))
        self.assertTrue(
            all(row["prediction_status"] == "ERROR_REF_MISMATCH" for row in mismatch)
        )
        self.assertTrue(all(row["delta_code_model"] == "ERROR" for row in mismatch))
        self.assertIn("WARNING: chr1:1001:C:G", standard_error)
        self.assertEqual(standard_error.count("WARNING: chr1:1001:C:G"), 1)

    def test_reference_requires_fai_index(self) -> None:
        fasta = self.root / "unindexed.fa"
        fasta.write_text(">chr1\nAAAA\n")
        with self.assertRaisesRegex(FileNotFoundError, "FASTA index not found"):
            vp.ReferenceValidator(fasta)

    def test_bundled_indel_table_has_all_exact_rows(self) -> None:
        data = vp.PredictionData(TOOL_DIR / "data")
        observed_rows = 0
        observed_keys = set()
        with gzip.open(TOOL_DIR / "data" / "MIC_indels.csv.gz", "rt", newline="") as handle:
            for row in csv.DictReader(handle):
                observed_rows += 1
                key = (row["variant"], row["event"], row["var_MIC_region"])
                self.assertNotIn(key, observed_keys)
                observed_keys.add(key)
                self.assertIn(row, data.indels_by_variant[row["variant"]])
        self.assertEqual(observed_rows, 21_966)
        self.assertEqual(len(observed_keys), 21_966)
        self.assertEqual(
            len({row["event"] for rows in data.indels_by_variant.values() for row in rows}),
            632,
        )
        self.assertEqual(len(data.events), 636)
        self.assertEqual(len(data.scale_by_event), 636)

    def test_known_notebook_snv_uses_st12_rescaling(self) -> None:
        self.assertEqual(
            vp.OUTPUT_COLUMNS[:13],
            [
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
            ],
        )
        vcf = self.write_vcf(
            "known.vcf", "chr1\t1051531\tknown\tC\tT\t.\tPASS\t.\n"
        )
        output_path = self.root / "known.tsv"
        standard_error = io.StringIO()
        with redirect_stderr(standard_error):
            exit_code = vp.main(
                [
                    "--input-vcf",
                    str(vcf),
                    "--output",
                    str(output_path),
                    "--data-dir",
                    str(TOOL_DIR / "data"),
                ]
            )
        self.assertEqual(exit_code, 0, standard_error.getvalue())
        with output_path.open(newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(len(rows), 1)
        row = rows[0]
        self.assertEqual(row["event"], "HsaEX0003157")
        self.assertEqual(row["var_MIC_region"], "c1")
        self.assertEqual(row["MIC_dist"], "-285")
        data = vp.PredictionData(TOOL_DIR / "data")
        expected = (
            vp.logit2(float(row["mt_code_score"]))
            - vp.logit2(float(row["wt_code_score"]))
        ) / data.scale_by_event[row["event"]]
        self.assertTrue(math.isclose(float(row["delta_code_model"]), expected, abs_tol=1e-14))
        self.assertFalse(
            math.isclose(float(row["delta_code_model"]), 0.0184489136621909, abs_tol=1e-8)
        )
        self.assertEqual(list(row), vp.OUTPUT_COLUMNS)

    def test_all_100_gnomad_alleles_match_expected_output(self) -> None:
        example = TOOL_DIR / "examples" / "gnomad_v4.1_microexon_100.vcf"
        expected = TOOL_DIR / "examples" / "gnomad_v4.1_microexon_100.expected.tsv"
        regions = set()
        events = set()
        novel_events = set()
        variant_types = []
        flank_variants = 0
        strongly_negative = 0
        with pysam.VariantFile(str(example)) as source:
            records = list(source)
        self.assertEqual(len(records), 100)
        for record in records:
            self.assertEqual(set(record.filter.keys()), {"PASS"})
            self.assertEqual(len(record.alts or ()), 1)
            self.assertIn("AF", record.info)
            self.assertIn("AC", record.info)
            self.assertIn("AN", record.info)
            regions.update(record.info["MIC_REGION"])
            events.update(record.info["MIC_EVENT"])
            if "new_mic" in record.info["MIC_GROUP"]:
                novel_events.update(record.info["MIC_EVENT"])
            if set(record.info["MIC_REGION"]) & {"up", "dn"}:
                flank_variants += 1
            if (
                record.info["MIN_DELTA_CODE_MODEL"] is not None
                and record.info["MIN_DELTA_CODE_MODEL"] <= -0.5
            ):
                strongly_negative += 1
            variant_types.append(vp.classify_variant(record.ref, record.alts[0]))
        self.assertEqual(regions, {"c1", "c2", "up", "dn", "ex"})
        self.assertEqual(len(events), 100)
        self.assertGreaterEqual(len(novel_events), 20)
        self.assertGreaterEqual(flank_variants, 95)
        self.assertGreaterEqual(strongly_negative, 80)
        self.assertEqual(variant_types.count("SNV"), 89)
        self.assertEqual(variant_types.count("Indel"), 11)

        output = self.root / "gnomad.tsv"
        standard_error = io.StringIO()
        with redirect_stderr(standard_error):
            exit_code = vp.main(
                [
                    "--input-vcf",
                    str(example),
                    "--output",
                    str(output),
                    "--data-dir",
                    str(TOOL_DIR / "data"),
                ]
            )
        self.assertEqual(exit_code, 0, standard_error.getvalue())
        self.assertEqual(output.read_bytes(), expected.read_bytes())
        self.assertIn("records=100", standard_error.getvalue())
        self.assertIn("ALT_alleles=100", standard_error.getvalue())
        self.assertIn("omitted_outside_coverage=0", standard_error.getvalue())
        with output.open(newline="") as handle:
            output_rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(list(output_rows[0]), vp.OUTPUT_COLUMNS)
        self.assertEqual(len({row["event"] for row in output_rows}), 100)
        self.assertGreaterEqual(
            sum(
                row["prediction_status"] == "FOUND"
                and float(row["delta_code_model"]) <= -0.5
                for row in output_rows
            ),
            80,
        )


if __name__ == "__main__":
    unittest.main()
