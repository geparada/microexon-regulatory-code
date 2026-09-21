# Microexon Code

## Installation

### Prerequisites

#### Hardware
- A x64 machine running Linux (tested on Ubuntu 18.04). Microexon Code does not run on macOS or Windows.
- An Nvidia GPU capable of CUDA 9.0+
- An internet connection as some required data is downloaded on demand
- If you want to run the code on large-scale data, you may need to run the code on a cluster or in the cloud. Snakemake is compatible with many schedulers and cloud APIs. See https://snakemake.readthedocs.io for details.

#### Software
Install Snakemake according to the instructions at https://snakemake.readthedocs.io/en/stable/getting_started/installation.html. This requires a Conda-based Python3 installation.

In short:

```shell
$ curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
$ bash Miniforge3-$(uname)-$(uname -m).sh
$ bash
```

Next, install Snakemake into a fresh environment and activate it:
```shell
$ mamba create -c conda-forge -c bioconda -n snakemake snakemake=8.14 
$ mamba activate snakemake
```

You will also need to add the `microexon_code/` directory to the local PYTHONPATH so the local packages can be found:
```shell
$ export PYTHONPATH=`pwd`  # or substitute the actual path to microexon_code/
```

## Overview
Microexon Code is a lightweight machine learning model using gradient boosted trees that is built on top of third-party feature detectors. To train the model or for inference, features need to be computed on the inputs. We use Snakemake to orchestrate this workflow. Some components require some inputs (e.g. model snapshots) which are supplied in the `resources/` folder.

The scripts supplied in this workflow can be used to train the model (even though we provide fully trained models that can be used for inference), replicate the analysis presented in the manuscript (such as _in silico_ saturation mutagenesis) and run inference on use supplied datasets.

## Usage
Microexon Code is run by executing Snakemake commands. Snakemake will then build a computational graph and run only those steps of the workflow required to obtain the desired output. The output can be specified by either supplying a named rule or supplying an _output_ file name.

### Inference for wild-type sequences or individual variants
In this mode, you may score sequences from a reference genome and optionally apply individual variants to each sequence. This is designed for a small number of variants which are applied to the wild-type sequence on the fly. Each input can have different variants applied. We use this mode to score SNVs such as our saturation _in silico_ mutagenesis experiments. 

To run inference on your own inputs you need to follow these steps exactly. The data format as well as the file names and paths must be exactly as specified as otherwise, Snakemake cannot deduce the correct workflow. For examples of correct inputs, see `resources/inputs/`

Create an input file in gzipped tab format with the following columns:
  - `event`: An event ID you can choose, e.g. `HsaEX6035820`
  - `chrom`: Chromosome in UCSC format (e.g. `chr5`, `chrX`)
  - `strand`: `+` or `-`
  - `upIntStart`, `upIntEnd`, `dnIntStart`, `dnIntEnd`: Genomic coordinates of the upstream and downstream introns of the microexon, 1-based indexing, e.g. `198639342       198692346       198692374       198696711`
  - `geneName`: Gene name, e.g. `PTPRC`
  - `lengthDiff`: The length of the microexon, e.g. `27`
  - `label`: Not used for inference, can be set to any value, e.g. `hi`
  - `variant` (optional column): Comma-delimited list of variants (e.g. `chr10:26770836:T:A`) to apply to the sequence. Variants must be non-overlapping. 

Examples of wild-type input tables are available under [`resources/inputs/wt/`](resources/inputs/wt/). For variant scoring, prepare your own table using the optional `variant` column described above; a separate variant-input example directory is not bundled.

Store your own input file under `resources/inputs/{task}/{prefix}.{species}.tab.gz`, where you can choose your own values for `{task}` and `{prefix}` and specify the desired genome build for `{species}`. Do not use any dots (`.`) in the prefix. An example of a valid path is `resources/inputs/wt/MicEvents.hg38.tab.gz`.
 
You can run inference on all input file, including the supplied examples, in `resources/inputs` via the command
```shell
  $ snakemake --cores {N} --use-conda --conda-not-block-search-path-envvars predictions
``` 
Substitute `{N}` for the number of cores on your machine or the number of parallel jobs you'd like to run.

Your predictions will be output in `results/predictions/{model_name}/{task}/predictions.{prefix}.{species}.csv.gz`

We will automatically output predictions for three different models which you can find under `{model_name}`

These are:
  - `known_microexons`: This model was trained on known microexons from seven mammalian species
  - `known_microexons_hg38`: This model was only trained on human known microexons
  - `novel_microexons`: This model was trained known microexons from seven mammalian species and the human novel microexons we detected in our study

### Saturation _in silico_ mutagenesis
You can run our saturation _in silico_ mutagenesis experiment with the command
```shell
  $ snakemake --cores {N} --use-conda --conda-not-block-search-path-envvars predict_saturation_in_silico_mutagenesis
```

This will generate a set of variants that mutate every nucleotide 150 nt upstream of a microexon, in a microexon, or 150 nt downstream of a microexon to every possible other nucleotide for all previously known human and mouse microexons.

The results will appear in `results/predictions/{model}/saturation_mutagenesis/predictions.{event}_variants.{species}.csv.gz`, where `{model}` is each of the three models above, `{event}` is all known microexons and `{species}` is `hg38` or `mm10`.

### QKI motif _in silico_ mutagenesis predictions
You can also run our _in silico_ mutagenesis experiment targeting QKI binding motifs. This script will find all QKI binding motifs of the pattern `ACTAA` or `ACTAAY` and mutate every nucleotide to every possible value. We do this for all known human and mouse microexons and generate predictions using each of our three models.

You can run this task with the command
```shell
  $ snakemake --cores {N} --use-conda --conda-not-block-search-path-envvars qki_predictions
```

Results will appear in `results/predictions/known_microexons/qki/predictions.qki_split_{split}_{motif}.{genome}.csv.gz` where `{split}` is a subset from 0 to 9.

### Inference for personal genomes
When we make predictions on microexon splicing for whole personal genomes such as our analysis of ASD-affected cohorts in the MSSNG and SSC datasets, we create a consensus genome file from the individuals VCF file and lift over the microexon coordinates to this new genome.

To run predictions with a personal genome, you need to place a bgzipped VCF file in `resources/personal_genomes/vcf/{dataset}/{sample}.vcf.gz` where you can choose `{dataset}` and `{sample}`. In our case `{dataset}` would be `MSSNG` or `SSC` and `{sample}` would correspond to the individual's sample id.

A sample-containing demonstration VCF is not bundled for this workflow. Supply your own authorized input VCF. For a runnable example of precomputed variant-score retrieval instead, see [`VariantPredictions/`](../VariantPredictions/); that example does not construct personal genomes.

The quality filter we apply to the VCF file is defined in [`workflow/rules/personal_genome_predictions.smk`](workflow/rules/personal_genome_predictions.smk#1):
```shell
FILTER="PASS" & 
FORMAT/DP>=10 & 
((GT="het" & (TYPE="snp" | TYPE="mnp") & GQ>=99 & AD[:1]/(FORMAT/DP)>=0.3 & AD[:1]/(FORMAT/DP)<0.8) | 
 (GT="het" & TYPE="indel" & GQ>=90 & AD[:1]/(FORMAT/DP)>=0.3 & AD[:1]/(FORMAT/DP)<0.8) | 
 (GT="AA" & GQ>=25 & AD[:1]/(FORMAT/DP)>=0.8))
```
Depending on your VCF file and the headers defined therein, you may need to modify this filter.

The personal-genome workflow expects its reference microexon table at `resources/inputs/personal_genomes/MicEvents_known_novel.hg38.tab.gz`. This file is not bundled; supply it in the input-table format described above before running the workflow.

After preparing the required inputs, run predictions on the VCF files under `resources/personal_genomes/vcf/{dataset}/{sample}.vcf.gz` with
```shell
  $ snakemake --cores {N} --use-conda --conda-not-block-search-path-envvars personal_genomes
```

The outputs will appear in ``results/predictions/{model_name}/personal_genomes/predictions.MicEvents_known_novel.[{dataset}--{sample}].csv.gz``.
