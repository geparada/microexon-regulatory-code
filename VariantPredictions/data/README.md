# VariantPredictions data

This directory contains the five gzip-compressed CSV files used by `../variant_predictions.py`. Together, they contain single-nucleotide substitution scores, wild-type scores, precomputed indel effects, event-level scaling factors, and genomic annotations for 636 microexon events. It also contains `microexon_metadata.hg38.csv`, a reference table describing each human microexon in this release and in Supplementary Table 2 of our paper. It is not used by the script, but can help users interpret its results. Researchers can use the files directly in custom pipelines by joining them through the shared `event` identifier. All genomic coordinates refer to GRCh38/hg38.

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

## `microexon_metadata.hg38.csv`

This table describes 1,858 human microexons: the 344 known microexons evaluated by the code model and the 1,514 human microexons in Supplementary Table 2 (ST2) of our paper (1,222 New and 292 New Neuronal). It records each microexon's annotation status, whether the code model evaluated and selected it, and whether it may need additional validation. Join it to the other tables using `event`. Microexon coordinates use the same 0-based, half-open convention as `microexon_events.hg38.csv.gz`.

| Column | Description |
|---|---|
| `microexon_id` | `chrom_strand_start_end`; the ST2 `ME-ID` for new microexons. |
| `chrom`, `strand`, `start`, `end` | Microexon coordinates (0-based, half-open). |
| `length` | Microexon length in nt. |
| `gene` | Gene name. |
| `event` | Event identifier used in the code release: VastDB `HsaEX…` for known microexons, `putative_microexon_…` for new ones. |
| `catalogue_class` | `known`: previously annotated microexon (VastDB). `new_neuronal`: new microexon with neuronal-differential inclusion (ST2 "New Neuronal"). `new`: new microexon detected but not neuronal-differential (ST2 "New"). |
| `in_variant_predictions` | `yes` for the 636 microexons with precomputed code-model scores in this directory (344 known, 292 New Neuronal). |
| `code_selected` | `yes`: one of the 369 microexons the code model recognises as carrying the neuronal regulatory code (ISM cluster 1, silhouette ≥ 0.1, RP11 clone-named genes removed, as in the Methods). These are the microexons used in the code and genetics analyses. `no`: evaluated but not selected. `not_evaluated`: not scored by the code model (all ST2 "New" microexons). |
| `ism_cluster`, `ism_silhouette` | ISM cluster and silhouette width behind `code_selected`. |
| `possible_alt_splice_site` | Whether the microexon could also correspond to an alternative 5' or 3' splice site of a flanking exon (see below): `exact` (identical sequence), `near` (1 to 2 mismatches), `sequence_match_only` (identical sequence, alternative junction not annotated) or `no`. |
| `alt_site_side` | `3' (acceptor)` or `5' (donor)`: the flanking exon with the alternative site. |
| `alt_site_mismatches` | Mismatches between the microexon and the bases at the alternative site. |
| `microexon_seq`, `alt_site_seq` | Microexon sequence and alternative-site bases, in transcript orientation. |

### Microexons that may correspond to alternative splice sites

For a small number of short microexons, the microexon sequence matches the bases between two alternative splice sites of a flanking exon, located the same number of nucleotides apart. In these cases, the mRNA that includes the microexon has the same sequence, or nearly the same sequence, as the mRNA in which the flanking exon is spliced at an alternative 5' or 3' splice site. RNA-seq reads, and RT-PCR products, may therefore not distinguish the two. These events are marked in `possible_alt_splice_site` and may require additional validation, for example genomic mutation of the microexon splice sites, before being interpreted as microexon inclusion. Nearly all code-selected microexons, which were used in the downstream genetic analyses of the paper, do not belong to this class.
