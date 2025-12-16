======
Output 
======

Output Structure
----------------

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


Output Examples
===============


