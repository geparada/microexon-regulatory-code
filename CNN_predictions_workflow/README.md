# CNN_predictions_workflow

Snakemake workflow used in this study to predict how genetic variants affect **microexon
splicing**, using two convolutional neural network (CNN) splicing models — **Pangolin** and
**DeltaSplice** — that we adapted specifically for microexons and applied to autism (ASD)
cohort variants.

This folder contains **analysis code only** (no data, inputs, or results — see *Data* below).

---

## What we used this workflow for

For every microexon in our catalogue, and for every candidate variant, the workflow asks a
focused question:

> *How much does this variant change the predicted usage / strength of this microexon's
> 3′ and 5′ splice sites?*

We ran it in three modes:

1. **Variant scoring** (`*_indels*.py`) — score real variants (from the SSC, SPARK and MSSNG
   ASD cohorts) at each microexon, producing a per-variant Δ splicing score.
2. **In-silico mutagenesis (ISM)** (`*_ISM.py`, `run_tissue_ISM.py`) — systematically mutate
   every position in and around each microexon to map the sequence determinants of splicing
   (e.g. regulatory hexamers).
3. **Sibling comparison** (`*_sib.py`, `get_sib_event_genomes*.py`) — build per-individual
   event sequences and compare affected probands with their unaffected siblings.

The variant-scoring outputs feed the downstream ASD variant analyses in
[`DataAnalysis/`](../DataAnalysis/). The separate [`VariantPredictions/data/`](../VariantPredictions/data/)
bundle contains code-model predictions, not the Pangolin or DeltaSplice outputs from this workflow.

## The concept behind the personalised Pangolin / DeltaSplice scripts

Pangolin and DeltaSplice are published, genome-wide variant-effect predictors. Out of the box
they scan whole genes and return a single, collapsed effect score. For this study we did **not**
use them that way — we re-purposed them into microexon-centred, tissue-resolved predictors.
The custom scripts here differ from the stock tools in four deliberate ways:

- **Microexon-centred, not genome-wide.** Instead of scanning a gene, each script places the
  model's receptive field directly on a microexon's **3′ (acceptor)** and **5′ (donor)** splice
  sites, taken from our microexon annotation, and reads the prediction at that exact site.

- **Per-tissue predictions, kept separate (Pangolin).** Pangolin has four tissue models
  (Heart, Liver, Brain, Testis). We load all four and read **each tissue's model at its own
  output channel**, reporting the four tissues independently — rather than collapsing them into
  a single max-across-tissues score as the stock tool does. DeltaSplice predictions use its
  reference-informed splice-site-usage (SSU) head, ensemble-averaged across the model set.

- **Correct indel handling.** After an insertion/deletion is applied, the splice-site coordinate
  shifts; the window is **re-centred on the post-variant position** (adjusting by the length
  change). Variants that span the splice site are evaluated over a small range of candidate
  positions and the strongest effect is kept.

- **WT vs mutant, reported as a delta.** Each script predicts the wild-type score and the
  mutant score at the microexon splice site and reports the difference (ΔPSI for Pangolin,
  ΔSSU for DeltaSplice) per event, variant, tissue and splice site.

These scripts **adapt**, and depend on, the original published models — we do not re-train them:

- **Pangolin** — Zeng & Li, *Genome Biology* 2022. Fork used: <https://github.com/geparada/Pangolin>
- **DeltaSplice** — Xu et al., *Genome Research* 2024. Fork used: <https://github.com/geparada/DeltaSplice>

The neural-network architectures and model weights belong to those works; the scripts in
`scripts/` are our microexon/ASD adaptations of them.

---

## Layout

```
Snakefile                     Orchestrates the rules (variant scoring + ISM)
scripts/
  # --- Pangolin (per-tissue P(splice) at microexon splice sites) ---
  pangolin_indels2.py         Variant scoring — the version run by the Snakefile
  pangolin_indels.py          Earlier variant-scoring version
  run_tissue_ISM.py           Per-tissue in-silico mutagenesis over microexon regions
  pangolin_sib.py             Proband vs unaffected-sibling scoring
  pred_ssu_ME.py              Splice-site-usage prediction helper
  # --- DeltaSplice (reference-informed SSU at microexon splice sites) ---
  deltasplice_indels4.py      Variant scoring — the version run by the Snakefile
  deltasplice_indels{,2,3}.py Earlier variant-scoring versions
  deltasplice_indels_cpu.py   CPU build of the above
  deltasplice_ISM.py          In-silico mutagenesis
  deltasplice_sib.py          Proband vs sibling scoring
  # --- shared ---
  get_sib_event_genomes{,_TOTAL}.py   Build per-individual event sequences for the sibling analysis
```

The active Snakefile rules run `pangolin_indels2.py` and `deltasplice_indels4.py` (variant
scoring) and `run_tissue_ISM.py` / `deltasplice_ISM.py` (ISM). The numbered siblings are earlier
iterations kept for provenance.

## Requirements

- [Snakemake](https://snakemake.readthedocs.io/)
- Two conda environments, referenced by name in the `Snakefile` (create them yourself; the
  original environment definitions were not part of this project):
  - `pangolin4` — Pangolin + PyTorch, NumPy, Biopython, pyfastx
  - `deltasplice` — DeltaSplice + PyTorch, NumPy, Biopython, pyfasta
- The **Pangolin** and **DeltaSplice** packages **and their model weights** installed in the
  respective environments (see the forks linked above).

## Data — not included

None of the following are in this repository; provide them locally and update the paths in the
`Snakefile` (and the model path near the top of `deltasplice_indels4.py`, currently a
cluster path):

| Input | Where to get it |
|---|---|
| Reference genome **hg38** (`hg38.fa`) | Public (UCSC / GENCODE / Ensembl) |
| Microexon annotation (`known_and_novel_microexons...hg38.csv`) | From this publication |
| Per-event variant lists (`input_indels*/…`) | Derived from the cohorts below |
| **SSC / SPARK / MSSNG** variants & metadata | **Controlled access** — SSC/SPARK via [SFARI Base](https://base.sfari.org/), MSSNG via [Autism Speaks MSSNG](https://research.mss.ng/). Subject to data-use agreements; **not redistributable** |

> ⚠️ The autism-cohort genotype/metadata are controlled-access human-subjects data and are
> intentionally excluded.

## Running

After creating the environments, installing the models, providing the inputs and editing the
paths in the `Snakefile`:

```bash
snakemake --use-conda --cores <N>
```
