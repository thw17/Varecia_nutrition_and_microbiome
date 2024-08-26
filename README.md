# Varecia_nutrition_and_microbiome

Code used in Beeby et al. Climate and nutrition drive gut microbiome variation in a fruit-specialist primate

This pipeline uses Snakemake (verision 7.28.1) to run and Conda/Bioconda to manage software environments (one for Snakemake to launch the pipeline and another to manage/launch QIIME2 for individual steps/commands--both are located in the ```envs``` subdirectory).

For more information about these tools, see:
Bioconda: https://bioconda.github.io/
Snakemake: https://snakemake.readthedocs.io/en/stable/

The Snakemake environment needs to be loaded before launching the pipeline, and once launched, Snakemake will manage the QIIME2 environment.