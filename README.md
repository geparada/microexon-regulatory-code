# A microexon regulatory code reveals autism-linked genetic variation


This repository contains the scripts and model predictions used in the analysis of autism-associated genetic variation and its impact on brain microexons, as detailed in our paper "A microexon regulatory code reveals autism-linked genetic variation".

## Description

Our study presents an integrative analysis of an expanded repertoire of neuronal differential microexons and the development of a machine learning model that delineates sequence elements critical for regulating microexon splicing. We reveal a landscape of predominantly rare human genetic variation that convergently impacts a subset of microexons in genes with known and unknown links to nervous system biology and neurodevelopmental disorders.

## Contents

- [`microexon-code/`](microexon-code/): This directory contains the computational workflow to perform predictions.
- [`DataAnalysis/`](DataAnalysis/): Analysis notebooks and input tables for microexon characterization, model evaluation and autism-associated genetic variation.
- [`VariantPredictions/`](VariantPredictions/): A standalone tool to retrieve precomputed code-model scores for variants in a GRCh38 VCF, with bundled prediction data and an example VCF.
- [`VariantPredictions/data/`](VariantPredictions/data/): This directory contains code-model predictions for all evaluated single-nucleotide substitutions and ASD-associated indels, together with event-level scaling factors and genomic annotations for the covered microexons.

## Precomputed variant scores

`VariantPredictions/` provides `variant_predictions.py`, a standalone script that takes a GRCh38 VCF and retrieves precomputed code-model scores for single-nucleotide variants and previously evaluated indels across 636 microexons and their surrounding regions. The directory includes the prediction tables, a test VCF containing 100 gnomAD variants, and its expected output. The bundled tables can also be used directly in custom pipelines. Instructions for installation and use, along with coverage and data descriptions, are provided within the directory.


## License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.

## Citation

If you use this code or the bundled predictions, please cite the associated paper:

> Parada, G. E.&#42;, Bretschneider, H.&#42;, Li, J. D., Fine, J. L., Dupas, S. J., Ellis, J. D., Braunschweig, U., Bonnal, S., Dalal, T., Engchuan, W., Zarrei, M., Wei, W., Zafar, N., Hemberg, M., Attisano, L., Ellis, J., Irimia, M., Wainberg, M., Trost, B., Scherer, S. W., Morris, Q. D.† & Blencowe, B. J.† A microexon regulatory code reveals autism-linked genetic variation. *Nature* (2026). [https://doi.org/10.1038/s41586-026-11119-w](https://doi.org/10.1038/s41586-026-11119-w)

The Zenodo code archive cited in the paper is:

> Parada, G. E.&#42;, Bretschneider, H.&#42;, Li, J. D., Fine, J. L., Dupas, S. J., Ellis, J. D., Braunschweig, U., Bonnal, S., Dalal, T., Engchuan, W., Zarrei, M., Wei, W., Zafar, N., Hemberg, M., Attisano, L., Ellis, J., Irimia, M., Wainberg, M., Trost, B., Scherer, S. W., Morris, Q. D.† & Blencowe, B. J.† Code for “A microexon regulatory code reveals autism-linked genetic variation”. *Zenodo* (2026). [https://doi.org/10.5281/zenodo.21940006](https://doi.org/10.5281/zenodo.21940006)
