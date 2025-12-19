======
Output 
======

Output Structure
================

::

	epigeneticbutton/
	├── config/			# Location for the main config file and recommended location for sample files and target files
	├── data/			# Location for test material and examples (e.g. zm_structural_RNAs.fa.gz)
	├── Help/			# Location for help files (e.g. Help_structural_RNAs_database_with_Rfam)
	├── profiles/
	│	├── sge/		# Config file to run snakemake on a cluster managed by SGE
	│	└── slurm/		# Config file to run snakemake on a cluster managed by SLURM
	├── workflow/
	│	├── envs/		# Conda environment file for depencies
	│	├── rules/		# Snakemake files with data type analysis rules
	│	├── scripts/		# R scripts for plots
	│	└── snakefile		# main snakefile
	├── genomes/			# Genome directories created upon run
	│	└── {ref_genome}/	# Reference genome directories with sequence, annotation and indexes
	└── results/			# Results directories created upon run
		├── combined/		# Combined analysis results
		│	├── bedfiles/	# Peak calling results
		│	├── chkpts/	# Empty checkpoint files used for pipeline logic. Deleting them will trigger rerunning the corresponding analysis
		│	├── logs/	# Log files
		│	├── matrix/	# Data matrices
		│	├── plots/	# Visualization plots
		│	└── reports/	# Analysis reports 
		└── <env>/	# Data type specific directories
			├── chkpts/	# Empty checkpoint files used for pipeline logic. Deleting them will trigger rerunning the corresponding analysis
			├── fastq/	# Processed FASTQ files
			├── logs/	# Log files
			├── mapped/	# Mapped reads (bam)
			├── plots/	# Data type specific plots
			├── reports/	# QC reports
			├── tracks/	# Track files (bigwigs)
			└── */		# data-specific directories (e.g. 'peaks' for ChIP, 'peaks' and 'motifs' for TF, 'DEG' for RNA, 'DMRs' and 'methylcall' for mC, 'clusters' for sRNA)

Data-specific Output
====================

Histone ChIP-seq
----------------

Output tree
+++++++++++

::

	ChIP/
	├── chkpts/	# Empty checkpoint files used for pipeline logic. Deleting them will trigger rerunning the corresponding analysis
	├── fastq/	# Processed FASTQ files
	├── logs/	# Log files
	├── mapped/	# Mapped reads (bam)
	├── peaks/	# Peak files (MACS2 output) for each replicate, pseudo-replicate and merged biological replicates and selected peaks (shared by merged and both pseudo-replicates).
	├── plots/	# Fingerprints (IP vs Input for each IP sample), IDR if at least two biological replicates
	├── reports/	# QC reports and summary of mapping statistics and peak statistics
	└── tracks/	# Track files (bigwigs); log2FC of IP/Input for each rep and merged if at least 2 biological replicates

Mapping statistics
++++++++++++++++++

- Data for each sample::

	results/ChIP/reports/summary_ChIP_<paired>_mapping_stats_ChIP__<line>__<tissue>__<sample_type>__<replicate>__<ref_genome>.txt

- Summary table:: 

	results/combined/reports/summary_mapping_stats_<analysis_name>_ChIP.txt

- Plot:: 

	results/combined/plots/mapping_stats_<analysis_name>_ChIP.pdf

- Example:

.. image:: images/mapping_stats_epicc_ChIP.png

(the actual output is in pdf format)

Peak statistics
+++++++++++++++

- Data for each sample:: 

	results/ChIP/reports/summary_ChIP_peak_stats_ChIP__<line>__<tissue>__<sample_type>__<ref_genome>.txt

- Summary table:: 

	results/combined/reports/summary_peak_stats_<analysis_name>_ChIP.txt

- Plot:: 

	results/combined/plots/peak_stats_<analysis_name>_ChIP.pdf

(see TF ChIP-seq for an example)

Fingerprints
++++++++++++

Performed with Deeptools.

- Plot for each biological replicate:: 

	results/ChIP/plots/Fingerprint__final__<data_type>__<line>__<tissue>__<sample_type>__<replicate>__<ref_genome>.png

(see TF ChIP-seq for an example) 

IDR
+++

Performed with IDR.

- Plot for pairs of biological replicate::

	results/ChIP/plots/idr_<paired>__<data_type>__<line>__<tissue>__<sample_type>__<replicate1>_vs_<replicate2>__<ref_genome>.<narrow|broad>Peak.png

(see TF ChIP-seq for an example)

Upset Plot
++++++++++

Perfomed with ComplexUpset.

- Table of combined peaks for all histone ChIP-seq samples in the analysis::

	results/combined/bedfiles/combined_peaks__ChIP__<analysis_name>__<ref_genome>.bed

- Table of combined peaks for all histone ChIP-seq samples in the analysis annotated based on the closest gene::

	results/combined/bedfiles/annotated__combined_peaks__ChIP__<analysis_name>__<ref_genome>.bed

- Upset plot::

	results/combined/plots/Upset_combined_peaks__ChIP__<analysis_name>__<ref_genome>.pdf

(see Combined Output for an example)

TF ChIP-seq
-----------

Output tree
+++++++++++

::

	TF/
	├── chkpts/	# Empty checkpoint files used for pipeline logic. Deleting them will trigger rerunning the corresponding analysis
	├── fastq/	# Processed FASTQ files
	├── logs/	# Log files
	├── mapped/	# Mapped reads (bam)
	├── motifs/	# Motifs analysis with the MEME suite, one folder per selected and idr peaks (and per replicates if so chosen in the config file)
	├── peaks/	# Peak files (MACS2 output) for each replicate, pseudo-replicate and merged biological replicates and selected peaks (shared by merged and both pseudo-replicates).
	├── plots/	# Fingerprints (IP vs Input for each IP sample), IDR if at least two biological replicates
	├── reports/	# QC reports and summary of mapping statistics and peak statistics
	└── tracks/	# Track files (bigwigs); log2FC of IP/Input for each rep and merged if at least 2 biological replicates

Mapping statistics
++++++++++++++++++

- Data for each sample:: 

	results/TF/reports/summary_TF_<paired>_mapping_stats_<data_type>__<line>__<tissue>__<sample_type>__<replicate>__<ref_genome>.txt

- Summary table:: 
	
	results/combined/reports/summary_mapping_stats_<analysis_name>_TF.txt

- Plot::
	
	results/combined/plots/mapping_stats_<analysis_name>_TF.pdf

(see histone ChIP-seq for an example) 

Peak statistics
+++++++++++++++

- Data for each sample::

	results/TF/reports/summary_TF_peak_stats_<dat_type>__<line>__<tissue>__<sample_type>__<ref_genome>.txt

- Summary table:: 

	results/combined/reports/summary_peak_stats_<analysis_name>_TF.txt

- Plot:: 

	results/combined/plots/peak_stats_<analysis_name>_TF.pdf

- Example:

.. image:: images/peak_stats_epicc_TF.png

(the actual output is in pdf format)

Fingerprints
++++++++++++

Performed with Deeptools.

- Plot for each biological replicate:: 

	results/ChIP/plots/Fingerprint__final__<data_type>__<line>__<tissue>__<sample_type>__<replicate>__<ref_genome>.png

- Example:

.. image:: images/Fingerprint__final__TF_SUVH1__Col0__suvh1.1__IP__Rep1__ColCEN.png

IDR
+++

Performed with IDR.

- Plot for pairs of biological replicate::

	results/ChIP/plots/idr_<paired>__<data_type>__<line>__<tissue>__<sample_type>__<replicate1>_vs_<replicate2>__<ref_genome>.<narrow|broad>Peak.png

- Example:

.. image:: images/idr_se__TF_SUVH1__Col0__suvh1.1__IP__Rep1_vs_Rep2__ColCEN_peaks.narrowPeak.png

Motifs
++++++

Performed with the MEME suite.

- Full output from selected peaks (and idr peaks if available) for each sample::

	results/TF/motifs/selected_peaks__<data_type>__<line>__<tissue>__<sample_type>__<ref_genome>/meme/

- which includes:: 

	results/TF/motifs/selected_peaks__<data_type>__<line>__<tissue>__<sample_type>__<ref_genome>/meme/meme_out/meme.html

- Example:

.. image:: images/meme.png

(the actual output is html format, and others)

Upset Plot
++++++++++

Perfomed with ComplexUpset.

- Table of combined peaks for all TF ChIP-seq samples in the analysis::

	results/combined/bedfiles/combined_peaks__TF__<analysis_name>__<ref_genome>.bed

- Table of combined peaks for all TF ChIP-seq samples in the analysis annotated based on the closest gene::

	results/combined/bedfiles/annotated__combined_peaks__TF__<analysis_name>__<ref_genome>.bed

- Upset plot::

	results/combined/plots/Upset_combined_peaks__TF__<analysis_name>__<ref_genome>.pdf

(see Combined Output for an example)

RNA-seq
-------

Output tree
+++++++++++

::

	RNA/
	├── chkpts/	# Empty checkpoint files used for pipeline logic. Deleting them will trigger rerunning the corresponding analysis
	├── DEG/ # Differential Expression Analysis results. Contains count tables, list of differential expression genes for all pairwise comparisons, gene expression tables and RData object for plotting gene expression (see `usage - plotting differential expression`).
	├── fastq/	# Processed FASTQ files
	├── GO/	# Gene Ontology Analysis results (optional). Contains GO terms enriched in sets of DEGs uniquely UP- or DOWN-regulated in each sample, and in additional GO analysis (see `usage - GO analysis`)
	├── logs/	# Log files
	├── mapped/	# Mapped reads (bam) (and STAR output files)
	├── plots/	# Expression and GO analysis (optional)
	├── reports/	# QC reports and summary of mapping statistics and peak statistics
	└── tracks/	# Track files (bigwigs); plus and minus strand (still in positive values) CPM for each replicate and merged all replicates per sample

Mapping statistics
++++++++++++++++++

- Data for each sample::

	results/RNA/reports/summary_RNA_<paired>_mapping_stats_<data_type>__<line>__<tissue>__<sample_type>__<replicate>__<ref_genome>.txt

- Summary table:: 
	
	results/combined/reports/summary_mapping_stats_<analysis_name>_RNA.txt

- Plot::
	
	results/combined/plots/mapping_stats_<analysis_name>_RNA.pdf

(see histone ChIP-seq for an example) 

Differential Expression analysis
++++++++++++++++++++++++++++++++

Counts from STAR; analysis performed with EdgeR.

- Count data for each RNAseq sample::

	results/RNA/DEG/counts__<data_type>__<line>__<tissue>__<sample_type>__<replicate>__<ref_genome>.tab

- Summary tables for all RNAseq samples used for the analysis:: 
	
	results/RNA/DEG/counts__<analysis_name>__<ref_genome>.txt # Count data output by STAR
	results/RNA/DEG/samples__<analysis_name>__<ref_genome>.txt # Table of samples information for edgeR analysis
	results/RNA/DEG/genes_rpkm__<analysis_name>__<ref_genome>.txt # Table of gene expression values for all genes in all samples in Reads per Kilobase Million (RPKM)

- Output tables of differentially expressed genes (DEG) for each pairwise comparison:: 
	
	results/RNA/DEG/FC_<analysis_name>__<ref_genome>__<line_sample1>__<tissue_sample1>_vs_<line_sample2>__<tissue_sample2>.txt # all genes in logFC sample1/sample2 and their differential statistics
	results/RNA/DEG/DEG_<analysis_name>__<ref_genome>__<line_sample1>__<tissue_sample1>_vs_<line_sample2>__<tissue_sample2>.txt # only DEGs

- Output summary tables of DEGs for all pairwise comparisons:: 

	results/RNA/DEG/summary_DEG_stats__<analysis_name>__<ref_genome>.txt # number of differential expressed genes in all pairwise comparisons and uniquely regulated in each sample
	results/RNA/DEG/unique_DEGs__<analysis_name>__<ref_genome>.txt # list of genes uniquely regulated in each sample

- Rdata object for plotting expression levels::

	results/RNA/DEG/ReadyToPlot__<analysis_name>__<ref_genome>.RData

- Global output from the differential analysis::

	results/combined/plots/BCV_RNAseq_<analysis_name>_<ref_genome>.pdf # Biological Coefficient of Variation of all genes
	results/combined/plots/MDS_RNAseq_<analysis_name>_<ref_genome>_d12.pdf # Multidimensional scaling of all the samples on the first two dimensions, with dots instead of labels
	results/combined/plots/MDS_RNAseq_<analysis_name>_<ref_genome>_d12_labs.pdf # Multidimensional scaling of all the samples on the first two dimensions, with labels instead of dots
	results/combined/plots/MDS_RNAseq_<analysis_name>_<ref_genome>_d23.pdf # Multidimensional scaling of all the samples on the first two dimensions, with dots instead of labels
	results/combined/plots/MDS_RNAseq_<analysis_name>_<ref_genome>_d23_labs.pdf # Multidimensional scaling of all the samples on the first two dimensions, with labels instead of dots

- Examples:

.. image:: images/BCV_RNAseq_epicc_ColCEN.png

.. image:: images/MDS2.png

- Heatmap of all DEGs across all samples::
	
	results/combined/plots/Heatmap_RNAseq_cpm__<analysis_name>__<ref_genome>.pdf # all gene expression normalized by count per million
	results/combined/plots/Heatmap_RNAseq_zscore__<analysis_name>__<ref_genome>.pdf # each gene normalized by Z-score

- Example:

.. image:: images/Heatmap_RNAseq_cpm__epicc__ColCEN.png

(the actual output is in pdf format)

- Plots of expression level in all samples for the top 100 DEGs (if present)::
	
	results/combined/plots/plot_expression__<analysis_name>__<ref_genome>__unique_DEGs.pdf

- Example:

.. image:: images/RNAseq_expression.png

Gene Ontology analysis
++++++++++++++++++++++

Performed with rrvgo and TopGO.

- List of Gene Ontology (GO) terms and corresponding Gene IDs (GIDs) enriched in the DEGs uniquely UP- and DOWN-regulated in each sample::

	results/RNA/GO/topGO_DOWN_in_<line>__<tissue>_BP_GOs.txt # Biological Process (BP) GO terms enriched in genes only DOWN-regulated in this sample 
	results/RNA/GO/topGO_DOWN_in_<line>__<tissue>_BP_GIDs.txt # genes in the Biological Process (BP) GO terms enriched in genes only DOWN-regulated in this sample 
	results/RNA/GO/topGO_DOWN_in_<line>__<tissue>_MF_GOs.txt # Molecular Function (MF) GO terms enriched in genes only DOWN-regulated in this sample
	results/RNA/GO/topGO_DOWN_in_<line>__<tissue>_MF_GIDs.txt # genes in the Molecular Function (MF) GO terms enriched in genes only DOWN-regulated in this sample
	results/RNA/GO/topGO_UP_in_<line>__<tissue>_BP_GOs.txt # Biological Process (BP) GO terms enriched in genes only UP-regulated in this sample
	results/RNA/GO/topGO_UP_in_<line>__<tissue>_BP_GIDs.txt # genes in the Biological Process (BP) GO terms enriched in genes only UP-regulated in this sample 
	results/RNA/GO/topGO_UP_in_<line>__<tissue>_MF_GOs.txt # Molecular Function (MF) GO terms enriched in genes only UP-regulated in this sample
	results/RNA/GO/topGO_UP_in_<line>__<tissue>_MF_GIDs.txt # genes in the Molecular Function (MF) GO terms enriched in genes only DOWN-regulated in this sample 
	
- Corresponding plots::

	results/RNA/plots/topGO_DOWN_in_<line>__<tissue>_BP_treemap.pdf # Treemap of simplified BP terms in DOWN-regulated genes in this sample
	results/RNA/plots/topGO_DOWN_in_<line>__<tissue>_MF_treemap.pdf # Treemap of simplified MF terms in DOWN-regulated genes in this sample
	results/RNA/plots/topGO_UP_in_<line>__<tissue>_BP_treemap.pdf # Treemap of simplified BP terms in UP-regulated genes in this sample
	results/RNA/plots/topGO_UP_in_<line>__<tissue>_MF_treemap.pdf # Treemap of simplified MF terms in UP-regulated genes in this sample

If not enough terms are enriched, these plots might not be created.

- Example (BP_DOWN):

.. image:: images/topGO_DOWN_in_Col0__suvh13_BP_treemap.png

(the actual output is in pdf format)

small RNA-seq
-------------

Output tree
+++++++++++

::

	sRNA/
	├── chkpts/	# Empty checkpoint files used for pipeline logic. Deleting them will trigger rerunning the corresponding analysis
	├── clusters/ # Clusters and differential analysis when all samples are analyzed together, on de novo identified clusters, on all genes and on all TEs (optional)
	├── fastq/	# Processed FASTQ files
	├── logs/	# Log files
	├── mapped/	# Subfolders of ShortStack output for each replicate
	├── reports/	# QC reports and summary of mapping statistics and peak statistics
	└── tracks/	# Track files (bigwigs); plus and minus strand (still in positive values) CPM for each replicate and merged all replicates per sample for each size chosen (default, 21, 22, 23 and 24nt)

Mapping statistics
++++++++++++++++++

- Data for each sample::

	results/sRNA/reports/sizes_stats__<data_type>__<line>__<tissue>__<sample_type>__<replicate>__<ref_genome>.txt

- Summary table:: 
	
	results/combined/reports/summary_sizes_stats_<analysis_name>_sRNA.txt

- Plot::
	
	results/combined/plots/srna_sizes_stats_<analysis_name>_sRNA.pdf # all sizes found in the raw data
	results/combined/plots/srna_sizes_stats_zoom_<analysis_name>_sRNA.pdf # zoom to chosen sizes (default 21 to 24nt)

- Example::

.. image:: images/srna_sizes_stats_epicc_sRNA.png

Cluster and Differential Expression analysis
++++++++++++++++++++++++++++++++++++++++++++

Counts from ShortStack; analysis performed with EdgeR.

- ShortStack analysis on each replicate::

	results/mapped/<data_type>__<line>__<tissue>__<sample_type>__<replicate>__<ref_genome>/ # output folder from ShortStack with all cluster results and alignement files
	results/mapped/<data_type>__<line>__<tissue>__<sample_type>__<replicate>__<ref_genome>/clusters.bed # simplified bed-file of clusters for downstream analyses

- For all samples in the analysis, two runs will be performed by default (three with the optional TE analysis), which will create the same output
- The full ShortStack output will be located in these folders::
	
	results/clusters/<analysis_name>__<ref_genome>__on_new_clusters/ # Identifying small RNA clusters, normalizing accross all samples
	results/clusters/<analysis_name>__<ref_genome>__on_all_genes/ # Limiting mapping to all genes, normalizing accross all samples
	results/clusters/<analysis_name>__<ref_genome>__on_all_TEs/ # Limiting mapping to all TEs, normalizing accross all samples (optional)

- Each folder (e.g. on new clusters) will also contain the files required for edgeR analysis and output from the differential analysis, following the same pattern than for DEG analysis of RNAseq::
	results/clusters/<analysis_name>__<ref_genome>__on_new_clusters/counts_for_edgeR.txt # Count data for edgeR analysis
	results/clusters/<analysis_name>__<ref_genome>__on_new_clusters/samples_for_edgeR.txt # Table of samples information for edgeR analysis
	results/clusters/<analysis_name>__<ref_genome>__on_new_clusters/FC_<line_sample1>__<tissue_sample1>_vs_<line_sample2>__<tissue_sample2>.txt # log Fold Change between each pairs of samples at all clusters
	results/clusters/<analysis_name>__<ref_genome>__on_new_clusters/DEG_<line_sample1>__<tissue_sample1>_vs_<line_sample2>__<tissue_sample2>.txt # only differentially regulated clusters between each pair of samples
	results/clusters/<analysis_name>__<ref_genome>__on_new_clusters/unique_DEGs.txt # list of clusters uniquely regulated in each sample

- For each differential analysis (e.g. on new clusters), similar output than for RNAseq DEGs will be generated, following the naming pattern ``sRNA_<analysis_name>_<ref_genome>__on_new_clusters`` pattern. It includes::
	
	results/sRNA/reports/summary_DEG_stats__<analysis_name>__<refgenome>__on_new_clusters.txt # number of differential regulated clusters in all pairwise comparisons and uniquely regulated in each sample
	results/combined/plots/BCV_RNAseq_<analysis_name>_<ref_genome>.pdf # Biological Coefficient of Variation of all genes
	results/combined/plots/MDS_RNAseq_<analysis_name>_<ref_genome>_<d12|d12_labs|d23|d23_labs>.pdf # Multidimensional scaling, all four versions
	results/combined/plots/Heatmap_sRNA_<cpm|zscore>__<analysis_name>__<ref_genome>__on_new_clusters.pdf # expression values accross all differentially regulated clusters by count per million and zscore.

(See RNAseq for examples)

Upset Plot
++++++++++

Perfomed with ComplexUpset.

- Table of combined clusters identified in at least one of the small RNA replicates in the analysis, split by chosen sizes::

	results/combined/bedfiles/combined_clusters__sRNA__<analysis_name>__<ref_genome>.bed

- Table of combined clusters annotated based on the closest gene::

	results/combined/bedfiles/annotated__combined_clusters__sRNA__<analysis_name>__<ref_genome>.bed

- Upset plot::

	results/combined/plots/Upset_combined_clusters__sRNA__<analysis_name>__<ref_genome>.pdf

- Example:

.. image:: images/Upset_combined_clusters__sRNA__epicc__ColCEN.png

DNA methylation
---------------

Combined Output
===============

Upset Plot
++++++++++

Perfomed with ComplexUpset.

- Table of combined peaks for all TF and histone ChIP-seq samples in the analysis::

	results/combined/bedfiles/combined_peaks__all_chip__<analysis_name>__<ref_genome>.bed

- Table of combined peaks for all TF and histone ChIP-seq samples in the analysis annotated based on the closest gene::

	results/combined/bedfiles/annotated__combined_peaks__all_chip__<analysis_name>__<ref_genome>.bed

- Upset plot::

	results/combined/plots/Upset_combined_peaks__all_chip__<analysis_name>__<ref_genome>.pdf

- Example:

.. image:: images/Upset_combined_peaks__all_chip__epicc__ColCEN.png

