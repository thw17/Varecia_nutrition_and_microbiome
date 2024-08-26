import os

# Runtime values
very_short = "6:00:00"
medium = "12:00:00"
day = "24:00:00"
long = "48:00:00"

metadata_file = "vv_microbiome_metadata_20221116_no_underscores.tsv"

rule all:
	input:
		"results/vv-table_250.qzv",
		"results/vv-rep-seqs_250.qzv",
		"results/vv-rep-seqs_250-aligned-rep-seqs.qza",
		"results/vv-rep-seqs_250_taxonomy.qzv",
		"results/vv-rep-seqs_250_taxonomy_barplot.qzv",
		"results/vv-rep-seqs_250_min2features_alphararefaction.qzv",
		# expand(
		# 	"diversity/results_alpha_depth{rarefaction}_vv-rep-seqs_250_min2features.txt", rarefaction=["1000", "3000", "6000", "10000"]),
		"features/results_featurefasta_vv-rep-seqs_250",
		"features/results_taxonomy_vv-rep-seqs_250/taxonomy.tsv",
		"features/results_featurecounts_vv-table_250/vv-table_250.tsv"

rule import_reads:
	input:
		"READS/"
	output:
		"results/vv-demux-paired-end.qza"
	params:
		threads = 4,
		mem = 16,
		t = medium
	conda:
		"envs/qiime2-2023.2-py38-linux-conda.yml"
	shell:
		"qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path {input} --input-format CasavaOneEightSingleLanePerSampleDirFmt --output-path {output}"

rule cutadapt:
	input:
		"results/vv-demux-paired-end.qza"
	output:
		"results/vv-demux-paired-end-trimmed.qza"
	params:
		threads = 4,
		mem = 16,
		t = medium
	conda:
		"envs/qiime2-2023.2-py38-linux-conda.yml"
	shell:
		"qiime cutadapt trim-paired --i-demultiplexed-sequences {input} --p-cores 4 --p-front-f GTGYCAGCMGCCGCGGTAA --o-trimmed-sequences {output} --verbose"

rule merge_pairs:
	input:
		"results/vv-demux-paired-end-trimmed.qza"
	output:
		"results/vv-demux-paired-end-trimmed-joined.qza"
	params:
		threads = 4,
		mem = 16,
		t = medium
	conda:
		"envs/qiime2-2023.2-py38-linux-conda.yml"
	shell:
		"qiime vsearch merge-pairs --i-demultiplexed-seqs {input} --o-merged-sequences {output}"

rule quality_filter:
	input:
		"results/vv-demux-paired-end-trimmed-joined.qza"
	output:
		seqs = "results/vv-demux-paired-end-trimmed-joined-filtered.qza",
		stats = "results/vv-demux-paired-end-trimmed-joined-filter-stats.qza"
	params:
		threads = 4,
		mem = 16,
		t = medium
	conda:
		"envs/qiime2-2023.2-py38-linux-conda.yml"
	shell:
		"qiime quality-filter q-score --i-demux {input} --o-filtered-sequences {output.seqs} --o-filter-stats {output.stats}"

rule deblur_denoise:
	input:
		"results/vv-demux-paired-end-trimmed-joined-filtered.qza"
	output:
		rep_seqs = "results/vv-rep-seqs_250.qza",
		table = "results/vv-table_250.qza",
		stats = "results/vv-deblur-stats_250.qza"
	params:
		threads = 4,
		mem = 16,
		t = medium
	conda:
		"envs/qiime2-2023.2-py38-linux-conda.yml"
	shell:
		"qiime deblur denoise-16S --i-demultiplexed-seqs {input} --p-trim-length 250 --p-sample-stats --o-representative-sequences {output.rep_seqs} --o-table {output.table} --o-stats {output.stats}"

rule summarize_initial_feature_table:
	input:
		"results/vv-table_250.qza"
	output:
		"results/vv-table_250.qzv"
	params:
		threads = 2,
		mem = 8,
		t = very_short
	conda:
		"envs/qiime2-2023.2-py38-linux-conda.yml"
	shell:
		"qiime feature-table summarize --i-table {input} --o-visualization {output}"

rule tabulate_initial_seqs:
	input:
		"results/vv-rep-seqs_250.qza"
	output:
		"results/vv-rep-seqs_250.qzv"
	params:
		threads = 2,
		mem = 8,
		t = very_short
	conda:
		"envs/qiime2-2023.2-py38-linux-conda.yml"
	shell:
		"qiime feature-table tabulate-seqs --i-data {input} --o-visualization {output}"

rule align_to_tree_mafft_fasttree:
	input:
		"results/vv-rep-seqs_250.qza"
	output:
		alignment = "results/vv-rep-seqs_250-aligned-rep-seqs.qza",
		masked_alignment = "results/vv-rep-seqs_250-masked-aligned-rep-seqs.qza",
		tree = "results/vv-rep-seqs_250-unrooted-tree.qza",
		rooted_tree = "results/vv-rep-seqs_250-rooted-tree.qza",
	params:
		threads = 4,
		mem = 16,
		t = medium
	conda:
		"envs/qiime2-2023.2-py38-linux-conda.yml"
	shell:
		"qiime phylogeny align-to-tree-mafft-fasttree --i-sequences {input} --o-alignment {output.alignment} --o-masked-alignment {output.masked_alignment} --o-tree {output.tree} --o-rooted-tree {output.rooted_tree}"

rule download_silva:
	output:
		"results/silva-138-99-nb-classifier.qza"
	params:
		web_address = "https://data.qiime2.org/2023.5/common/silva-138-99-nb-classifier.qza",
		threads = 1,
		mem = 4,
		t = very_short
	conda:
		"envs/qiime2-2023.2-py38-linux-conda.yml"
	shell:
		"wget {params.web_address} -O {output}"

rule feature_classifier:
	input:
		silva = "results/silva-138-99-nb-classifier.qza",
		seqs = "results/vv-rep-seqs_250.qza"
	output:
		"results/vv-rep-seqs_250_taxonomy.qza"
	params:
		threads = 4,
		mem = 16,
		t = medium
	conda:
		"envs/qiime2-2023.2-py38-linux-conda.yml"
	shell:
		"qiime feature-classifier classify-sklearn --i-classifier {input.silva} --i-reads {input.seqs} --o-classification {output}"

rule tabulate_classified_seqs:
	input:
		"results/vv-rep-seqs_250_taxonomy.qza"
	output:
		"results/vv-rep-seqs_250_taxonomy.qzv"
	params:
		threads = 2,
		mem = 8,
		t = very_short
	conda:
		"envs/qiime2-2023.2-py38-linux-conda.yml"
	shell:
		"qiime metadata tabulate --m-input-file {input} --o-visualization {output}"

rule barplot_taxonomy:
	input:
		table = "results/vv-table_250.qza",
		tax = "results/vv-rep-seqs_250_taxonomy.qza"
	output:
		"results/vv-rep-seqs_250_taxonomy_barplot.qzv"
	params:
		metadata = metadata_file,
		threads = 2,
		mem = 8,
		t = very_short
	conda:
		"envs/qiime2-2023.2-py38-linux-conda.yml"
	shell:
		"qiime taxa barplot --i-table {input.table} --i-taxonomy {input.tax} --m-metadata-file {params.metadata} --o-visualization {output}"

rule remove_singletons:
	input:
		"results/vv-table_250.qza"
	output:
		"results/vv-rep-seqs_250_min2features.qza"
	params:
		threads = 2,
		mem = 8,
		t = very_short
	conda:
		"envs/qiime2-2023.2-py38-linux-conda.yml"
	shell:
		"qiime feature-table filter-features --i-table {input} --p-min-samples 2 --o-filtered-table {output}"

rule summarize_nosingleton_feature_table:
	input:
		"results/vv-rep-seqs_250_min2features.qza"
	output:
		"results/vv-rep-seqs_250_min2features.qvz"
	params:
		threads = 2,
		mem = 8,
		t = very_short
	conda:
		"envs/qiime2-2023.2-py38-linux-conda.yml"
	shell:
		"qiime feature-table summarize --i-table {input} --o-visualization {output}"

rule alpha_rarefaction:
	input:
		table = "results/vv-rep-seqs_250_min2features.qza",
		phylo = "results/vv-rep-seqs_250-rooted-tree.qza"
	output:
		"results/vv-rep-seqs_250_min2features_alphararefaction.qzv"
	params:
		threads = 2,
		mem = 8,
		t = very_short,
		m = metadata_file
	conda:
		"envs/qiime2-2023.2-py38-linux-conda.yml"
	shell:
		"qiime diversity alpha-rarefaction --i-table {input.table} --i-phylogeny {input.phylo} --m-metadata-file {params.m} --p-max-depth 20000 --p-min-depth 1 --p-steps 20 --o-visualization {output}"

rule alpha_diversity_analyses:
	input:
		table = "results/vv-rep-seqs_250_min2features.qza",
		phylo = "results/vv-rep-seqs_250-rooted-tree.qza"
	output:
		shannon = "diversity/min2features_{rarefaction}_core_metrics_results/shannon_vector.qza",
		faith = "diversity/min2features_{rarefaction}_core_metrics_results/faith_pd_vector.qza"
	params:
		threads = 2,
		mem = 8,
		t = very_short,
		m = metadata_file,
		depth = "{rarefaction}",
		o = "diversity/min2features_{rarefaction}_core_metrics_results"
	conda:
		"envs/qiime2-2023.2-py38-linux-conda.yml"
	shell:
		"qiime diversity core-metrics-phylogenetic --i-table {input.table} --i-phylogeny {input.phylo} "
		"--m-metadata-file {params.m} --p-sampling-depth {params.depth} "
		"--o-shannon-vector {output.shannon} --o-faith-pd-vector {output.faith} "
		"--o-rarefied-table {params.o}/table.qza --o-observed-features-vector {params.o}/observed_features.gza "
		"--o-evenness-vector {params.o}/evenness.qza --o-unweighted-unifrac-distance-matrix {params.o}/unweighted_unifrac.qza "
		"--o-weighted-unifrac-distance-matrix {params.o}/weighted_unifrac.qza --o-jaccard-distance-matrix {params.o}/jaccard.qza "
		"--o-bray-curtis-distance-matrix {params.o}/bray-curtis.qza "
		"--o-unweighted-unifrac-pcoa-results {params.o}/unweighted_unifrac_pcoa.qza --o-weighted-unifrac-pcoa-results {params.o}/weighted_unifrac_pcoa.qza "
		"--o-jaccard-pcoa-results {params.o}/jaccard_pcoa.qza "
		"--o-bray-curtis-pcoa-results {params.o}/bray-curtis_pcoa.qza --o-unweighted-unifrac-emperor {params.o}/unweighted_unifrac.qzv "
		"--o-weighted-unifrac-emperor {params.o}/weighted_unifrac.qzv "
		"--o-jaccard-emperor {params.o}/jaccard.qzv --o-bray-curtis-emperor {params.o}/bray-curtis.qzv"

rule metadata_merge_alpha:
	input:
		shannon = "diversity/min2features_{rarefaction}_core_metrics_results/shannon_vector.qza",
		faith = "diversity/min2features_{rarefaction}_core_metrics_results/faith_pd_vector.qza"
	output:
		"diversity/results_alpha_depth{rarefaction}_vv-rep-seqs_250_min2features.qzv"
	params:
		threads = 2,
		mem = 8,
		t = very_short,
		m = metadata_file
	conda:
		"envs/qiime2-2023.2-py38-linux-conda.yml"
	shell:
		"qiime metadata tabulate --m-input-file {params.m} --m-input-file {input.shannon} --m-input-file {input.faith} --o-visualization {output}"

rule export_merged_alpha:
	input:
		"diversity/results_alpha_depth{rarefaction}_vv-rep-seqs_250_min2features.qzv"
	output:
		directory("diversity/results_alpha_depth{rarefaction}_vv-rep-seqs_250_min2features")
	params:
		threads = 2,
		mem = 8,
		t = very_short
	conda:
		"envs/qiime2-2023.2-py38-linux-conda.yml"
	shell:
		"qiime tools export --input-path {input} --output-path {output}"

rule export_feature_counts:
	input:
		"results/vv-rep-seqs_250_min2features.qvz"
	output:
		directory("features/results_featurecounts_vv-rep-seqs_250_min2features")
	params:
		threads = 2,
		mem = 8,
		t = very_short
	conda:
		"envs/qiime2-2023.2-py38-linux-conda.yml"
	shell:
		"qiime tools export --input-path {input} --output-path {output}"

rule export_output_feature_fasta:
	input:
		"results/vv-rep-seqs_250.qza"
	output:
		directory("features/results_featurefasta_vv-rep-seqs_250")
	params:
		threads = 2,
		mem = 8,
		t = very_short
	conda:
		"envs/qiime2-2023.2-py38-linux-conda.yml"
	shell:
		"qiime tools export --input-path {input} --output-path {output}"

rule export_output_feature_table:
	input:
		"results/vv-table_250.qza"
	output:
		"features/results_featurecounts_vv-table_250/feature-table.biom"
	params:
		threads = 2,
		mem = 8,
		t = very_short
	conda:
		"envs/qiime2-2023.2-py38-linux-conda.yml"
	shell:
		"qiime tools export --input-path {input} --output-path {output}"

rule export_taxonomy:
	input:
		"results/vv-rep-seqs_250_taxonomy.qza"
	output:
		"features/results_taxonomy_vv-rep-seqs_250/taxonomy.tsv"
	params:
		threads = 2,
		mem = 8,
		t = very_short
	conda:
		"envs/qiime2-2023.2-py38-linux-conda.yml"
	shell:
		"qiime tools export --input-path {input} --output-path {output}"

rule biom_convert:
	input:
		"features/results_featurecounts_vv-table_250/feature-table.biom"
	output:
		"features/results_featurecounts_vv-table_250/vv-table_250.tsv"
	params:
		threads = 2,
		mem = 8,
		t = very_short
	conda:
		"envs/qiime2-2023.2-py38-linux-conda.yml"
	shell:
		"biom convert -i {input} -o {output} --to-tsv"