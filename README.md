# Varecia_nutrition_and_microbiome

Code used in Beeby et al. Climate and nutrition drive gut microbiome variation in a fruit-specialist primate

This pipeline uses Snakemake (verision 7.28.1) to run and Conda/Bioconda to manage software environments (one for Snakemake to launch the pipeline and another to manage/launch QIIME2 for individual steps/commands--both are located in the ```envs``` subdirectory).

For more information about these tools, see:
Bioconda: https://bioconda.github.io/
Snakemake: https://snakemake.readthedocs.io/en/stable/

The Snakemake environment needs to be loaded before launching the pipeline, and once launched, Snakemake will manage the QIIME2 environment.

Note that the pipeline is expecting the fastq files to be in a subdirectory called "READS". This can be updated in the ```import_reads``` rule's input directive (should be line 27). The pipeline also expects the fastq files to be named as the were coming off of the sequencer. Any updates/changes to read names need to be reflected in the metadata file's first column.