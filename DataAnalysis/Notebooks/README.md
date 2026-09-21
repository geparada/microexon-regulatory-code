# Data-analysis notebooks

This directory contains analysis notebooks for **A microexon regulatory code reveals autism-linked genetic variation** by Parada, Bretschneider and colleagues.

## Analysis topics

Each topic index groups notebooks by figure or analysis and lists figure and statistics reproductions and figure workflows before supporting and historical work.

- [Microexon discovery and characterization](discovery-and-characterization/)
- [NMD and developmental regulation](nmd-and-development/)
- [Code model and regulatory sequences](code-model-and-regulatory-sequences/)
- [Experimental and external validation](experimental-and-external-validation/)
- [Autism-associated variant analyses](autism-variant-analyses/)
- [Gene function and protein context](gene-function-and-protein-context/)
- [Supplementary tables and public resources](supplementary-tables-and-public-resources/)

## Notebook roles

- `figures/`: figure-analysis workflows.
- `reproductions/`: focused notebooks that reconstruct specific panels or statistics.
- `supporting/`: input preparation, methods and additional statistics.
- `tables/`: supplementary-table preparation.
- `archive/`: earlier analyses retained for historical context.

A notebook has one home even when it contributes to several analyses. Start with the relevant topic index; historical notebooks can use earlier data versions or analysis settings.

## Running the notebooks

See the [parent directory](../) for environment setup. These are research analysis notebooks, not self-contained demonstrations. Analyses using controlled-access cohort data require independently authorized access to those inputs.

The topic folders organize the notebooks; they are not new analysis working directories. Before executing a notebook, set its kernel working directory to a **writable local analysis workspace** containing the expected inputs and output directories. For example, run `%cd /path/to/your/analysis-workspace` in a scratch cell. Resolve the notebook's relative paths, including `Results/`, `CONTROLLED_ACCESS/`, `EXTERNAL_DATA/` and any `../` paths, from that workspace. Adapt absolute paths to your installation. Merely opening a notebook from its new folder does not configure these paths.

Notebook code and saved results have been retained during the directory reorganization. For a self-contained variant-retrieval demonstration, use [VariantPredictions/](../../VariantPredictions/).

## Related resources

- [microexon-code/](../../microexon-code/): model training, inference and in silico mutagenesis.
- [CNN_predictions_workflow/](../../CNN_predictions_workflow/): complementary variant-scoring workflows.
- [VariantPredictions/data/](../../VariantPredictions/data/): substitution and indel predictions, scaling factors and microexon coordinates.
- [Input/](../Input/): bundled analysis inputs.
