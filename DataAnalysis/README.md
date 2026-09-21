# Data Analysis for Microexon Code Paper

This directory contains analysis scripts for "A microexon regulatory code reveals autism-linked genetic variation".


## Notebooks

The [notebook index](Notebooks/) links to figure analyses and related resources. Some notebooks have cleared outputs; others retain reviewed outputs. Participant identifiers that were redacted are represented by placeholders such as `REDACTED_PARTICIPANT`.

Paths beginning with `CONTROLLED_ACCESS/` refer to inputs or derived files that are not distributed here. Running those analyses requires independently authorized access to the underlying data and adjustment of the paths for your local setup. `EXTERNAL_DATA/` marks other external resources that are not bundled. The original analysis logic is retained; these notebooks are not self-contained demonstrations.

Before sharing an executed notebook, review its outputs for individual-level information. From the repository root, run `python tools/scan_cohort_ids.py` as an additional check. It checks known identifier patterns and saved state in cohort-related notebooks, including downloaded Git LFS files and decompressed tables. It does not replace a data-sharing review or inspect text within images.


## Setup

To run these notebooks, you can create a virtual environment that contains most of the necessary dependencies using the following command:

```sh
mamba create --name DataAnalysis numpy pandas pybedtools pybigwig pysam r-data.table r-dbplyr r-ggplot2 r-ggsignif r-reshape rpy2 scipy seaborn bedtools bioconductor-biobase bioconductor-annotationdbi bioconductor-iranges bioconductor-pcamethods r-cowplot
```
Then activate environment 

```sh
conda activate DataAnalysis
```

To open and run these notebook files, you can use [JupyterLab](https://jupyter.org/), which provides an interactive environment for running and editing Jupyter notebooks.

An exhaustive list of all the libraries included in the virtual environment we used to carry out these analyses and their corresponding versions is found in the `env/environment_clone.yml` file.

##  Input Tables
Some input tables necessary to run these notebooks are available in the `Input/` directory. Controlled-access participant data are not included and must be obtained through the relevant data provider under the applicable access conditions.
