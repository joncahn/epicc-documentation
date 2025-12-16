Usage
=====

Running the pipeline
--------------------

1. To run the pipeline locally:

``snakemake --use-conda --conda-frontend conda --cores 12``

2. To run the pipeline on a HPC-slurm (using sbatch):

``snakemake --profile profiles/slurm``

If you do not want all the snakemake output (very talkative), instead of using `--quiet` I would recommmend redirecting it to a log and putting the run in the background:

``snakemake --profile profiles/slurm > epigeneticbutton.log 2>&1 &``

3. Other option: To run the pipeline on a UGE cluster (using qsub):

``mkdir hpclogs``

``snakemake --profile profiles/uge``

*The commands for the clusters are specific to the CSHL environement. If using a profile, make sure these parameters are adapted to your cluster too or edit accordingly. See operating systems help.*

4. Optional: to test the pipeline, consider generating a DAG first to make sure your samplefiles and parameters work:

``snakemake --dag | dot -Tpng > dag.png``

*Even if snakemake is launched on a cluster with a profile option, the run will output a lot on the terminal. It is recommended to launch the command from a screen, to start it from a script submitted to the cluster, or to put the command in the background (which will still output snakemake commands but allows further action).*

*For full understanding of snakemake capabilities and options: https://snakemake.readthedocs.io/en/stable/*

Intermediate Target Rules
--------------------------

Snakemake allows intermediate files or rules to be targeted instead of the whole pipeline. If a specific output file is specificied after the snakemake run command, only the steps leading to this file will be performed, not the rest of the pipeline.
For example, to only generate the plot of mapping statistics for the ChIP samples of the analysis name "mysamples" (and all the intermediate files required for it):
``snakemake --cores 1 results/combined/plots/mapping_stats_mysamples_ChIP.pdf``

In addition to the additional output plots below, two rules can be specified as targets for intermediates outputs:

	- ``map_only``: Only performs the alignement of all samples. It returns bam files, QC files and mapping metrics. Usage: ``snakemake --cores 1 map_only``

	- ``coverage_chip``: Creates bigwig files of coverage for all ChIP samples. The binsize is by default 1bp (can be updated in the config file ``chip_tracks: binsize: 1``). Usage: ``snakemake --cores 1 coverage_chip``

Additional Output Options
--------------------------

Below is a list of *cool* outputs that can be generated once whole pipeline ran once. You'll find a basic structure for how to tell snakemake to generate them, feel free to replace the --cores 1 with the HPC profile you would rather use.

**1. Plotting RNAseq expression levels on target genes**
++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Given a list of genes (and optional labels), it will plot the expression levels in all the different samples in the samplefile and analysis name defined. Genes uniquely differentially regulated in one sample versus one or more samples are color coded. It is based on a Rdata file created during the Differential Expression analysis.
To run it, edit the config file with the target gene list file (``rnaseq_target_file``: 1 column list of genes ID that must match the gtf file of the reference genome used, optional second column for gene labels, additional columns can be present but will not be used) and a corresponding label (``rnaseq_target_file_label``: name which will be included in the name of the output pdf) and run the following command, replacing <analysis_name>, <ref_genome> and <rnaseq_target_file_label> with wanted values:

::

  snakemake --cores 1 results/RNA/plots/plot_expression__<analysis_name>__<ref_genome>__<rnaseq_target_file_label>.pdf

Note that the separators between the variables are two underscores next to each other "__", except in "plot_expression".
An example where <analysis_name>="test_smk" and <ref_genome>="TAIR10", while setting the target file and its label "my_genes_of_interests" directly in the snakemake command:

::

  snakemake --cores 1 results/RNA/plots/plot_expression__test_smk__TAIR10__my_genes_of_interests.pdf --config rnaseq_target_file="data/target_genes.txt" rnaseq_target_file_label="my_genes_of_interests"

Output is a single pdf file named ``results/RNA/plots/plot_expression__<analysis_name>__<ref_genome>__<rnaseq_target_file_label>.pdf`` where each gene of the list is on an individual page.

**2. Performing GO analysis on target genes**
+++++++++++++++++++++++++++++++++++++++++++++

Given a file containing a list of genes to do GO analysis on, and optionally a background file (default to all genes in the reference genome), it will perform Gene Ontology analysis.
By default, GO is not performed since it requires manual input to build a database. To activate it, ``GO`` needs to be switched to ``true`` in the config file, and the files to make the GO database should be defined in the config file ``gaf_file`` and ``gene_info_file`` below the corresponding reference genome. See `Help_GO <https://github.com/joncahn/epigeneticbutton/blob/main/Help/Help_Gene_Ontology>`__ for more details on how to create the GO database.
To run it, edit the config file with the target gene list file (``rnaseq_target_file``: 1 column list of genes ID that must match the gtf file of the reference genome used, optional second column for gene labels, additional columns can be present but will not be used) and a corresponding label (``rnaseq_target_file_label``: name which will be included in the name of the output files) and run the following command, replacing <analysis_name>, <ref_genome> and <rnaseq_target_file_label> with wanted values:

::

  snakemake --cores 1 results/RNA/GO/TopGO__<analysis_name>__<ref_genome>__<rnaseq_target_file_label>.done

Note that the separators between the variables are two underscores next to each other "__".
An example where <analysis_name>="test_smk" and <ref_genome>="ColCEN", while setting the target file and its label "my_genes_of_interests" directly in the snakemake command:

::

  snakemake --cores 1 results/RNA/GO/TopGO__test_smk__ColCEN__my_genes_of_interests.done --config rnaseq_target_file="data/target_genes.txt" rnaseq_target_file_label="my_genes_of_interests"

Output are two pdf files, one for the biological process terms ``results/RNA/plots/topGO_<rnaseq_target_file_label>_BP_treemap.pdf`` and one for the molecular function terms ``results/RNA/plots/topGO_<rnaseq_target_file_label>_MF_treemap.pdf``. Corresponding tables listing the terms enriched for each gene of the ``rnaseq_target_file`` are also generated at ``results/RNA/GO/topGO_<rnaseq_target_file_label>_<BP|MF>_<GOs|GIDs>.txt`` for a focus on the GO terms or the GIDs, respectively.

**3. Finding motifs on target regions**
+++++++++++++++++++++++++++++++++++++++

Given a bed file containing different regions, it will perform a motifs analysis with meme.
By default motifs analysis is only performed on the final selected TF peak files (``motifs: true`` in the config file). Edit to ``allrep: true`` in the config file for motifs analysis to be performed on all replicates and pairwise idr peaks if available. A plant motifs database is used by default for tomtom. Download the appropriate file from JASPAR and replace its name in the config file ``jaspar_db`` and change the ``motifs_ref_genome`` to match the samples.
To run the analysis:

::

  snakemake --cores 1 results/TF/chkpts/motifs__<motif_target_file_label>.done

Note that the separator is two underscores next to each other "__".
An example running the pipeline on a slurm hpc, for regions from <ref_genome>="ColCEN", while setting the target file and its label "my_genes_of_interests" directly in the snakemake command:

::

  snakemake --profile profiles/slurm results/TF/chkpts/motifs__my_regions_of_interests.done --config motifs_target_file="data/target_peaks.txt" motifs_target_file_label="my_regions_of_interests" motifs_ref_genome="ColCEN"

Output is the folder ``results/TF/<motif_target_file_label>`` containing a subdirectory called ``meme`` and potentially one called ``tomtom`` with all the results, as described in https://meme-suite.org/meme/index.html.
When setting ``motif_ref_genome``, it is safer to use a reference genome that has already been used in a run. Otherwise, it will be treated like the *ref_genome* of a sample, creating a fasta file in the ``genomes/<ref_genome>`` directory if a fasta file is found at ``ref_path``.
For the target file chosen ``motif_target_file``, if the regions are over 500bp, only the middle 400bp will be used.

**4. Performing sRNA differential analysis on regions**
+++++++++++++++++++++++++++++++++++++++++++++++++++++++

Given a bed or gff file, it will perform the small RNA analysis with shortstack followed by differential analysis with edgeR, using all the samples from the sample file but limiting the mapping and counts to the loci in the target file. Edit ``srna_target_file`` and ``srna_target_file_label`` in the config file. To run the analysis: 

::

  snakemake --cores 1 results/sRNA/clusters/<analysis_name>__<ref_genome>__on_<srna_target_file_label>/Counts.txt

Note that the separators between the variables are two underscores next to each other "__" except between "on" and "<srna_target_file_label>" where it's only one "_".
An example running the pipeline on a slurm hpc, <analysis_name>="test_smk" and <ref_genome>="ColCEN", while setting the target file and its label "miRNAs" directly in the snakemake command:

::

  snakemake --profile profiles/slurm results/sRNA/clusters/test_smk__ColCEN__on_miRNAs/Counts.txt --config sRNA_target_file="data/miRNA.gff" sRNA_target_file_label="miRNAs"

Output is the results folder from Shortstack limited to this loci file, followed by the differential cluster analysis with edgeR.

If you only want the results of Shortstack and not the differential analysis, limit the run to the rule ``analyze_all_srna_samples_on_target_file`` instead, by targeting: ``results/sRNA/clusters/<analysis_name>__<ref_genome>__on_<srna_target_file_label>/Counts.txt``.

The bed or gff file of regions **MUST HAVE** a header with a column called "Name" (the 4th column of a bed file or the 9th column of a gff3).

**5. Plotting heatmap on regions**
++++++++++++++++++++++++++++++++++

Given a bed file, it will plot a heatmap using deeptools.
Edit ``heatmap_target_file`` and ``heatmap_target_file_label`` in the config file. To run the analysis: 

::

  snakemake --cores 1 results/combined/plots/Heatmap__<matrix_param>__<env>__<analysis_name>__<ref_genome>__<target_name>.pdf

- the <matrix_param> can be ``regions`` for scaled regions, ``tss`` for reference point on the TSS or ``tes`` for reference point on the TES.
- the <env> correspond to the data types to include. Since mC requires different parameters, it has to be done independently. If you have several different data types including mC, and want the order of the regions to be maintained in the mC heatmap based on the other samples, use:

::

  snakemake --cores 1 results/combined/plots/Heatmap_sorted__<matrix_param>__mC__<analysis_name>__<ref_genome>__<target_name>.pdf

This will generate the heatmap for all the other samples first. If you want the regions sorted based on the mC samples only, use:

::

  snakemake --cores 1 results/combined/plots/Heatmap__<matrix_param>__mC__<analysis_name>__<ref_genome>__<target_name>.pdf

To make a heatmap will all the samples (excluding mC), use <env>="most". If you want to include mC samples (will probbaly *not work*) use <env>="all".

An example running the pipeline on a slurm hpc, with <analysis_name>="test_smk", <ref_genome>="ColCEN", <matrix_param>="regions", on all samples but mC <env>="most", while setting the target file and its label "interesting_genes" directly in the snakemake command:

::

  snakemake --profile profiles/slurm results/combined/plots/Heatmap__regions__most__test_smk__ColCEN__interesting_genes.pdf --config heatmap_target_file="data/target_genes.bed" heatmap_target_file_label="interesting_genes"

Output is a pdf file, or two if sorted heatmap for mC samples was generated.
By default, the heatmaps will be scaled by type (i.e. each ChIP mark, each TF, RNAseq, each sRNAseq size and each mC context on an appropriate scale based on the values in the heatmap). It can be changed to "default", where a single scale is used for the whole heatmap, or to "sample" where each sample is scaled individually. This can be changed in the config file ``heatmaps_scales``.
By default, the heatmaps are sorted based on "mean" of all samples accross all regions. This can be changed in the config file ``heatmaps_sort_options`` to "median" or to "none", keeping the regions in the order of the bedfile.
If the given bedfile is stranded, the heatmap will be done by splitting the regions into plus and minus strand for proper stranded data (RNAseq and sRNAs) values. If this is not the wanted behavior, disable ``stranded_heatmaps`` in the config file.
The color scheme of the heatmap is "seismic" for all samples and "Oranges" for mC. This can be changed manually in the config file ``heatmaps_plot_params``.
The size of the scaled regions ``middle`` (-m in deeptools), the size of the surrounding regions ``before`` (-b in deeptools) and ``after`` (-a in deeptools) and the binsize ``binsize`` (-bs in deeptools) can be edited in the config file in ``heatmaps`` for each <matrix_params>.

**6. Plotting metaplot profiles on regions**
++++++++++++++++++++++++++++++++++++++++++++

Given a bed file, it will plot a metaplot profile using deeptools.
Edit ``heatmap_target_file`` and ``heatmap_target_file_label`` in the config file. To run the analysis: 

::

  snakemake --cores 1 results/combined/plots/Profile__<matrix_param>__<env>__<analysis_name>__<ref_genome>__<target_name>.pdf

Similar to heatmap above for the <matrix_param> options.
Use <env>="all" to include all samples (mC and others).
Output is two pdf files, where the samples are grouped by regions or not.
By default, the heatmaps will be scaled by type (i.e. each ChIP mark, each TF, RNAseq, each sRNAseq size and each mC context on their appropriate scale based on the values in the heatmap). It can be changed to "default", where a single scale is used for the whole heatmap, or to "sample" where each sample is scaled individually. This can be changed in the config file ``heatmaps_scales``.
By default, the profiles represent the "mean" accross all regions. This can be changed in the config file ``profile_scale`` to "median".
By default, the type of plots are "lines". See deeptools documentation for other options and edit ``profiles_plot_params`` in the config file.
The size of the scaled regions ``middle`` (-m in deeptools), the size of the surrounding regions ``before`` (-b in deeptools) and ``after`` (-a in deeptools) and the binsize ``binsize`` (-bs in deeptools) can be edited in the config file in ``heatmaps`` for each <matrix_params>.

**7. Plotting browser screenshots on regions**
++++++++++++++++++++++++++++++++++++++++++++++

Given a region file, it will plot a browser screenshot using R packages.
Edit ``browser_target_file`` and ``browser_target_file_label`` in the config file. To run the analysis: 

::

  snakemake --cores 1 results/combined/plots/Browser_<target_name>__<env>__<analysis_name>__<ref_genome>.pdf

The target file is a bed-like file, with the following columns: 
+-----+-------+------+-------+---------+----------------------------+----------------------------+
| Chr | Start | End  | ID    | Binsize | Higlight_starts (optional) | Higlight_widths (optional) |
+=====+=======+======+=======+=========+============================+============================+
| chr1 | 1000 | 5000 | Peak1 | 1       | 3000,4000                  | 50,200                     |
+-----+-------+------+-------+---------+----------------------------+----------------------------+

Each region will be printed individually, and merged into a final PDF.
Hightlights columns are optional, and correspond to regions of the browser that will be highlighted for this specific region (boxed). As many highlights can be used in a comma-separated lists, the first highlight will be in blue and all the others in red. On the above example, the region to plot is chr1:1000-5000, using col6=3000,4000 and col7=50,200 will make a blue box higlighting chr1:3000-3050 and a red one highlighting chr1:4000:4200.
Use <env>="all" to include all samples, "most" for all data-types except mC, or any single environment for data type-specific browsers [all, most, ChIP, TF, RNA, sRNA, mC].
By default, no TE file is used. If you want to add TE annotations, supply a bed-file in the config file ``browser_TE_file``.

### **8. Rerunning a specific analysis**
++++++++++++++++++++++++++++++++++++++++

To rerun a specific analysis, force snakemake to recreate the target file, adding to the snakemake command: ``<target_file> --force``
e.g ``snakemake --cores 1 results/combined/plots/srna_sizes_stats_test_snakemake_sRNA.pdf --force``
If only the combined analysis is to be performed, and not everything else, delete all the chkpts files in ``results/combined/chkpts/`` as well as in the chkpt of each relevant environment ``results/<env>/chkpts/<env>_analysis__<analysis_name>__<ref_genome>.done``.
Changing parameters in the config file should trigger a rerun of the impacted samples.
Several instances of the epicc pipeline can thus be chained to change some parameters in a script, for example:

::

  for scale in default sample type; do
    snakemake --profile profiles/slurm results/combined/plots/Heatmap__regions__most__test_smk_${scale}__ColCEN__interesting_genes.pdf --config heatmap_target_file="data/target_genes.bed" heatmap_target_file_label="interesting_genes" heatmaps_scales=${scale}
  done

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


