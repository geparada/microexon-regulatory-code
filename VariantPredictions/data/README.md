# VariantPredictions data

This directory contains the five gzip-compressed CSV files used by `../variant_predictions.py`. Together, they contain single-nucleotide substitution scores, wild-type scores, precomputed indel effects, event-level scaling factors, and genomic annotations for 636 microexon events. Researchers can use the files directly in custom pipelines by joining them through the shared `event` identifier. All genomic coordinates refer to GRCh38/hg38.

## `MIC_ISM.gz`

This table contains 3,459,549 single-nucleotide substitutions evaluated by in silico mutagenesis (ISM) across the five model regions. Each row records the prediction for one substitution in one microexon event and region.

| Column | Description |
|---|---|
| `event` | Microexon event identifier used to join this table to the wild-type score and event metadata tables. |
| `region` | Evaluated model region: microexon (`ex`), upstream intronic region (`up`), downstream intronic region (`dn`), or a window around one of the two flanking splice sites (`c1` or `c2`). Upstream and downstream are defined relative to transcript orientation. |
| `rel_pos` | Position of the substituted nucleotide within the source ISM region coordinate system. Its value can differ from the `MIC_dist` convention used in Supplementary Table 12 of our paper. |
| `ref` | Reference nucleotide at the evaluated position. |
| `alt` | Alternate nucleotide introduced during in silico mutagenesis. |
| `mt_score` | Raw code-model probability after introducing the alternate nucleotide. |
| `delta_logit_score` | Difference between the alternate and wild-type scores after conversion to log2 odds, before event-level scaling. |
| `transformed_delta_scores` | Effect score from the source ISM analysis, calculated by scaling `delta_logit_score` within each event. Retained for reproducibility; `variant_predictions.py` calculates its output `delta_code_model` from the raw scores and the scaling factors reported with the paper. |

## `MIC_ISM_wt.gz`

This table contains one wild-type code-model score for each of the 636 events. Join it to `MIC_ISM.gz` using `event`.

| Column | Description |
|---|---|
| `event` | Microexon event identifier. |
| `wt_score` | Raw code-model probability for the unmodified, wild-type event sequence. |

## `MIC_indels.csv.gz`

Derived from the data reported in Supplementary Table 12 of our paper, this table contains 21,966 predictions for combinations of indels, microexon events, and regions, representing 20,935 distinct variant strings across 632 events. These indels were identified in ASD cohorts from the Simons Simplex Collection (SSC), SPARK, and MSSNG. An indel appears in multiple rows when it was evaluated for more than one event or region.

| Column | Description |
|---|---|
| `variant` | Exact indel representation in `chrom:position:REF:ALT` format, using a 1-based VCF position. |
| `event` | Microexon event identifier. |
| `var_MIC_region` | Evaluated region matched by the indel: `c1`, `c2`, `up`, `dn`, or `ex`. |
| `MIC_dist` | Position relative to the relevant microexon or splice-site anchor, using the convention in Supplementary Table 12. |
| `mt_code_score` | Raw code-model probability for the sequence carrying the indel. |
| `wt_code_score` | Raw code-model probability for the corresponding wild-type event sequence. |
| `logit_scaling_factor` | Event-specific factor used for the paper's main analyses, based on the original analysis variant set. |
| `delta_code_model` | Change in log2 odds divided by the original `logit_scaling_factor` used in the paper's main analyses. |
| `logit_scaling_factor_all_vars` | Event-specific factor recalculated using all raw code-model predictions in this table. |
| `delta_code_model_all_vars` | Change in log2 odds divided by `logit_scaling_factor_all_vars`, with values between -1 and 1. Returned by default by `variant_predictions.py`. |

The two delta columns use different scaling factors for each microexon event. `delta_code_model` uses the original factors from the paper's main analyses. Because this table includes additional variants, applying those factors can produce values below -1 or above 1. `delta_code_model_all_vars` uses factors recalculated from all raw code-model predictions in this table, so its values remain between -1 and 1.

For custom lookups, match the complete `variant` string and keep all matching rows, since an indel may have predictions for more than one microexon event or region.

## `MIC_event_scaling.csv.gz`

This table contains the two event-level scaling factors used by `variant_predictions.py` to calculate SNV effects. It has one row for each of the 636 microexon events considered by the data bundle.

| Column | Description |
|---|---|
| `event` | Microexon event identifier. |
| `logit_scaling_factor` | Event-specific factor used for the paper's main analyses. Selected with `--delta-scaling main-analysis`. |
| `logit_scaling_factor_all_vars` | Event-specific factor recalculated across the expanded variant set. Used by default. |

## `microexon_events.hg38.csv.gz`

This table defines the genomic coordinates and annotations of the 636 covered microexon events. Intron intervals use 0-based, half-open coordinates: the start is included and the end is excluded. Upstream and downstream are defined relative to transcript orientation; `Start` and `End` always denote the lower and upper genomic boundaries, respectively.

| Column | Description |
|---|---|
| `event` | Microexon event identifier used across all five bundled tables. |
| `chrom` | GRCh38 chromosome with a `chr` prefix. |
| `strand` | Transcriptional strand, `+` or `-`. |
| `upIntStart` | Lower genomic boundary of the transcript-upstream intron. |
| `upIntEnd` | Upper genomic boundary of the transcript-upstream intron. |
| `dnIntStart` | Lower genomic boundary of the transcript-downstream intron. |
| `dnIntEnd` | Upper genomic boundary of the transcript-downstream intron. |
| `geneName` | Gene symbol associated with the microexon. |
| `lengthDiff` | Microexon length in nucleotides. |
| `group` | Source annotation group: `known_mic` or `novel_mic`. In the output, `variant_predictions.py` uses `new_mic` for novel microexons to match the naming in Supplementary Table 12. |

When mapping ISM substitutions to genomic positions in a custom pipeline, combine these event coordinates with `region` and `rel_pos`, accounting for strand and the difference between 0-based interval coordinates and 1-based VCF positions. `variant_predictions.py` implements this conversion for VCF queries.
