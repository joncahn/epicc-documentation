=============
Configuration
=============

For new users, it is recommended to use the configuration `epicc-builder app <https://epicc-builder.streamlit.app/>`__ to validate your sample metadata file and choose analysis options.

Metadata samplefile
===================

Summary
-------

Prepare your sample metadata file (default to ``config/all_samples.tsv``) with the required 9 columns below (see below for more details specific to each data-type):

+-------------+--------+----------+---------------+-------------+------------+--------------+----------+--------------+
| *data_type* | *line* | *tissue* | *sample_type* | *replicate* | *seq_id*   | *fastq_path* | *paired* | *ref_genome* |
+=============+========+==========+===============+=============+============+==============+==========+==============+

  - Col1: *data_type*
      Type of data. Only takes one of these options: [RNAseq | ChIP | TF | mC | sRNA]
	  See below for how to further specify ChIP and TF samples.

  - Col2: *line*
      Sample line (e.g. ``B73``). It can also be used to define different genotypes (e.g. ``WT``, ``ddm1``, ``dnmt3a``), or other characteristics that vary between samples, such as collection time points.
  
  - Col3: *tissue* 
      Tissue type (e.g. ``Leaf``). It can also be used to define different genotypes (e.g.``WT``, ``ddm1``, ``dnmt3a``), or other characteristics that vary between samples, such as collection time points.
  
  - Col4: *sample_type* 
      Sample identifier. It depends on the type of data; for example for ChIP-seq data, it is used to differentiate the IP from the corresponding Input. See `Columns unique per data type` below for the different options.
  
  - Col5: *replicate*
      Replicate ID (e.g ``Rep1``, ``RepA``, ``plate56``). Different replicates must have the same *line* and *tissue* values to be analyzed together as replicates.
  
  - Col6: *seq_id*
      Sequencing ID. Unique identifier to identify the raw data. Use the corresponding SRR####### if downloading from SRA, or a unique string which is unique to         your sample in the fastq path below (but shared between Read1 and Read2 for paired-end data).
  
  - Col7: *fastq_path* 
      Path to FASTQ files. if downloading from SRA, use ``SRA``.
  
  - Col8: *paired*
      Whether the raw data is paired-end (``PE``) or single-end (``SE``). Only takes one of these two options: [PE | SE]
  
  - Col9: *ref_genome* 
      Reference genome name. Will need an entry in the configuration file.

A template can be found on the `epicc-builder app <https://epicc-builder.streamlit.app/>`__ and you can use it to confirm that your entries follow the epxected patterns.

Example
-------

A test example is below (the header is only indicative, and should not be present on the actual file), using data from `Cahn et al. 2024 <https://pubmed.ncbi.nlm.nih.gov/39632087/>`__ and `Lee et al. 2020 <https://pubmed.ncbi.nlm.nih.gov/32303559/>`__.

+-------------+--------+--------------+---------------+-------------+-------------+--------------+----------+--------------+
| *data_type* | *line* | *tissue*     | *sample_type* | *replicate* | *seq_id*    | *fastq_path* | *paired* | *ref_genome* |
+=============+========+==============+===============+=============+=============+==============+==========+==============+
ChIP          | Col0   | WT           | H3K27ac       | Rep1        | SRR27821885 | SRA          | PE       | ColCEN       |
+-------------+--------+--------------+---------------+-------------+-------------+--------------+----------+--------------+
ChIP          | Col0   | WT           | H3K4me1       | Rep1        | SRR27821852 | SRA          | PE       | ColCEN       |
+-------------+--------+--------------+---------------+-------------+-------------+--------------+----------+--------------+
ChIP          | Col0   | WT           | Input         | Rep1        | SRR27821932 | SRA          | PE       | ColCEN       |
+-------------+--------+--------------+---------------+-------------+-------------+--------------+----------+--------------+
mC            | Col0   | suvh1        | WGBS          | Rep1        | SRR27821959 | SRA          | SE       | ColCEN       |
+-------------+--------+--------------+---------------+-------------+-------------+--------------+----------+--------------+
mC            | Col0   | WT           | WGBS          | Rep1        | SRR27821907 | SRA          | SE       | ColCEN       |
+-------------+--------+--------------+---------------+-------------+-------------+--------------+----------+--------------+
mC            | Col0   | WT           | WGBS          | Rep2        | SRR27821907 | SRA          | SE       | ColCEN       |
+-------------+--------+--------------+---------------+-------------+-------------+--------------+----------+--------------+
RNAseq        | Col0   | suvh13       | RNAseq        | Rep1        | SRR27821840 | SRA          | PE       | ColCEN       |
+-------------+--------+--------------+---------------+-------------+-------------+--------------+----------+--------------+
RNAseq        | Col0   | suvh13       | RNAseq        | Rep2        | SRR27821839 | SRA          | PE       | ColCEN       |
+-------------+--------+--------------+---------------+-------------+-------------+--------------+----------+--------------+
RNAseq        | Col0   | suvh13       | RNAseq        | Rep3        | SRR27821838 | SRA          | PE       | ColCEN       |
+-------------+--------+--------------+---------------+-------------+-------------+--------------+----------+--------------+
RNAseq        | Col0   | WT           | RNAseq        | Rep1        | SRR27821967 | SRA          | PE       | ColCEN       |
+-------------+--------+--------------+---------------+-------------+-------------+--------------+----------+--------------+
RNAseq        | Col0   | WT           | RNAseq        | Rep2        | SRR27821966 | SRA          | PE       | ColCEN       |
+-------------+--------+--------------+---------------+-------------+-------------+--------------+----------+--------------+
RNAseq        | Col0   | WT           | RNAseq        | Rep3        | SRR27821918 | SRA          | PE       | ColCEN       |
+-------------+--------+--------------+---------------+-------------+-------------+--------------+----------+--------------+
TF_SUVH1      | Col0   | suvh1.1      | Input         | Rep1        | SRR27821931 | SRA          | SE       | ColCEN       |
+-------------+--------+--------------+---------------+-------------+-------------+--------------+----------+--------------+
TF_SUVH1      | Col0   | suvh1.1      | Input         | Rep2        | SRR27821934 | SRA          | SE       | ColCEN       |
+-------------+--------+--------------+---------------+-------------+-------------+--------------+----------+--------------+
TF_SUVH1      | Col0   | suvh1.1      | IP            | Rep1        | SRR27821933 | SRA          | SE       | ColCEN       |
+-------------+--------+--------------+---------------+-------------+-------------+--------------+----------+--------------+
TF_SUVH1      | Col0   | suvh1.1      | IP            | Rep2        | SRR27821935 | SRA          | SE       | ColCEN       |
+-------------+--------+--------------+---------------+-------------+-------------+--------------+----------+--------------+
TF_SUVH3      | Col0   | suvh3        | Input         | RepA        | SRR27821929 | SRA          | SE       | ColCEN       |
+-------------+--------+--------------+---------------+-------------+-------------+--------------+----------+--------------+
TF_SUVH3      | Col0   | suvh3        | IP            | RepA        | SRR27821930 | SRA          | SE       | ColCEN       |
+-------------+--------+--------------+---------------+-------------+-------------+--------------+----------+--------------+
sRNA          | Col0   | VLP_ddm1     | shRNA         | Rep1        | SRR8792540  | SRA          | SE       | ColCEN       |
+-------------+--------+--------------+---------------+-------------+-------------+--------------+----------+--------------+
sRNA          | Col0   | VLP_ddm1     | shRNA         | Rep2        | SRR8792538  | SRA          | SE       | ColCEN       |
+-------------+--------+--------------+---------------+-------------+-------------+--------------+----------+--------------+
sRNA          | Col0   | VLP_WT_inflo | shRNA         | Rep1        | SRR8792539  | SRA          | SE       | ColCEN       |
+-------------+--------+--------------+---------------+-------------+-------------+--------------+----------+--------------+

Columns common to all types of samples
--------------------------------------

- Col2: *line*
   Can be any information you want, such as ``Col0`` or ``WT`` to annotate and label samples

- Col3: *tissue* 
   Can be any information you want, such as ``leaf`` or ``mutant`` or ``6h_stress`` to annotate and label samples.
   The combination line x tissue will be the base for all comparisons (e.g ``WT_leaf`` vs ``WT_roots`` or ``Col0_control`` vs ``Ler_stress``)

- Col5: *replicate* 
   Any value to match the different replicates (e.g ``Rep1``, ``RepA``, ``1``). All the different replicates are merged for samples with the same *line (Col2)*, *tissue (Col3)* and *sample_type (Col4)*.

- Col6: *seq_id*
   Unique identifier to identify the raw data. If the data is deposited in SRA, it can be an SRR number (e.g. ``SRR27821931``) or a comma-delimited list of SRR numbers (e.g. ``SRR27821931,SRR27821932,SRR27821933``) without spaces if multiple fastq files should be merged into 1 biological replicate. If the data is local, it must be a unique identifier of the file in this folder (e.g. ``wt_k27``). This identifier should be shared by both 'R1' and 'R2' fastq files for paired-end data.

- Col7: *fastq_path*
   Either ``SRA`` if raw data to be downloaded from SRA (the SRR number should be used as *seq-id*), or the path to the directory containing the fastq file (e.g. ``/archive/fastq``), in which case the *seq_id* should be a unique identifier of the corresponding file in this folder (e.g. ``/archive/fastq/raw.reads.wt_k27.fastq.gz``)

- Col8: *paired*
   ``PE`` for paired-end data or ``SE`` for single-end data. PE samples should have two fastq files 'R1' and 'R2' at the location defined in * fastq_path (Col7)*, sharing the same *seq-id (Col6)* (e.g. ``/archive/fastq/raw.reads.wt_k27_R1.fastq.gz`` and ``/archive/fastq/raw.reads.wt_k27_R2.fastq.gz``)

- Col9: *ref_genome*
	Name of the reference genome to use for mapping (e.g ``tair10``). 
  	For each reference genome, a corresponding fasta, gff and gtf files are required. It can be a full path (including the extension) or relative to the main repo folder. These files can be gzipped. For example, if your sample file has ``B73_NAM`` as a *ref_genome (Col9)*, there must be this entry in the config file: 

::

	B73_NAM:
		fasta_file: path/to/B73.fasta	# can be .fa(.gz) or .fasta(.gz)
		gff_file: B73.gff	# can be .gff*(.gz)
		gtf_file: B73.gtf	# can be .gtf(.gz)

Other files specific to each reference genome are optional.
The GTF file can be created from a GFF file with cufflinks ``gffread -T <gff_file> -o <gtf_file>`` and check that ``transcript_id`` and ``gene_id`` are correctly assigned in the 9th column. The GFF file should have ``gene`` and ``exon`` in the 3rd column. All files can be gzipped (.gz extension).

Columns specific to each data type
----------------------------------

Histone ChIP-seq
^^^^^^^^^^^^^^^^^^^

- Col1: *data_type*	
	``ChIP`` or ``ChIP_<id>`` where ``<id>`` is an identifier to relate an IP sample to its corresponding input. Only necessary in case there are different inputs to be used for different IP samples that otherwise share the same *line (Col1)* and *tissue (Col2)* values.
   	For example: If you have H3K27ac IP samples which you want normalized to an H3 sample, and H4K16ac to be normalized to an H4 sample. Both H3 and H4 samples should be labeled ``Input`` in *sample_type (Col4*), so to differentiate them, use ``ChIP_H3`` and ``ChIP_H4`` for the *data_type* of these Inputs and for the H3K27ac and H4K16ac IPs, respectively.

Example:

+--------+------+----+---------+------+------------+----------+----+--------+
|ChIP_H3 | Col0 | WT | H3K27ac | Rep1 | wt_k27     | ./fastq/ | PE | ColCEN |
+--------+------+----+---------+------+------------+----------+----+--------+
|ChIP_H3 | Col0 | WT | IP      | Rep1 | wt_h3_ctrl | ./fastq/ | PE | ColCEN |
+--------+------+----+---------+------+------------+----------+----+--------+
|ChIP_H4 | Col0 | WT | H3K27ac | Rep1 | wt_h4k16   | ./fastq/ | PE | ColCEN |
+--------+------+----+---------+------+------------+----------+----+--------+
|ChIP_H4 | Col0 | WT | IP      | Rep1 | wt_h4_ctrl | ./fastq/ | PE | ColCEN |
+--------+------+----+---------+------+------------+----------+----+--------+

- Col4: *sample_type*
    Either ``Input`` to be used as a control (even if it is actually H3 or IgG pull-down), or the histone mark IP (e.g. ``H3K9me2``). If the mark is not already listed in the config file ``chip_callpeaks: peaktype:``, add it to the desired category (either narrow or broad peaks).
   
- Option:
      Differential nucleosome sensitivity (DNS-seq) can be analyzed with ``ChIP`` *data_type*, using ``MNase`` for the light digest and ``Input`` for the heavy digest.

Transcription factor ChIP-seq
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- Col1: *data_type*
	``TF_<tf_name>`` where ``<tf_name>`` is the name of the transcirption factor (e.g. for ``TB1`` data, use ``TF_TB1``). This name should be identical for the IP and its input, and for all replicates. Multiple TFs can be analyzed in parallel, each having its own set of IP and Input samples e.g. ``TF_<name1>`` and ``TF_<name2>``.
   
- Col4: *sample_type*
	Either ``Input`` or ``IP``. This works for transcription factors with narrow peaks (default). Use ``IPb`` for broad peaks.

RNA-seq
^^^^^^^^^^

- Col1: *data_type*
	``RNAseq``. No other options.

- Col4: *sample_type*
	``RNAseq``. No other options.

small RNA-seq
^^^^^^^^^^^^^^^^

- Col1: *data_type*
	``sRNA``. No other options.

- Col4: *sample_type*
	``sRNA``. Can also be ``smallRNA`` or ``shRNA``. Does not change the analysis but it is used in file names.

Whole Genome DNA methylation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- Col1: *data_type*
	``mC``. No other options.

- Col4: *sample_type*
	``mC``. Can also be ``WGBS``, ``ONT``, ``Pico`` or ``EMseq``. Not yet relevant except for the file names, but will be used in future release to define the type of data.

Configuration file
==================

Summary
-------

Update ``config/config.yaml`` with your paths and parameters.
  
- Sample file
	This is the full path to the file detailed above which contain your samples metadata. 
  
- Reference genome files 
	For each reference genome in the sample file (*ref_genome (Col9)*), enter the corresponding species, the full paths to a fasta file, a gene gff file, and a gene gtf file. Optionally, a path to a gaf file and gene info file can be given for Gene Ontology (GO) analysis (see `Extra outputs` and `Help GO <https://github.com/joncahn/epigeneticbutton/blob/main/Help/Help_Gene_Ontology>`__ for more details), an annotation file (bed format) for transposable elements for TE-centered analysis, and a fasta file of structural RNAs to be depleted from small RNA data (see `Help Rfam <https://github.com/joncahn/epigeneticbutton/blob/main/Help/Help_structural_RNAs_database_with_Rfam>`__ for more details).
example: 

::

	ColCEN:
		species: "thaliana"
		fasta_file: path/to/ColCEN.fasta	# can be .fa(.gz) or .fasta(.gz)
		gff_file: path/to/ColCEN.gff	# can be .gff*(.gz)
		gtf_file: path/to/ColCEN.gtf	# can be .gtf(.gz)
		gaf_file: "data/ColCEN_infoGO.tab.gz" # optional. can be gzipped or not.
  		gene_info_file: "data/ColCEN_genes_info.tab.gz" # optional. can be gzipped or not.
  		te_file: "path/to/ColCEN_TEs.bed.gz" # optional. can be gzipped or not.
  		structural_rna_fafile: "path/to/ColCEN_structural_RNAs.fa.gz" # optional.

- Analysis parameters / options:
	Works with default parameter, but more details can be found below, on the `epicc-builder app <https://epicc-builder.streamlit.app/>`__, or directly commented on the ``config/config.yaml`` file for more customization.
  
- Species-specific parameters:
	For each species in the reference genomes used, the number to used for the STAR index and the size of the genome are required.
	The NCBI ID, genus and go_database are only required for gene ontology (GO) analysis, if the option has been selected (See `Help GO <https://github.com/joncahn/epigeneticbutton/blob/main/Help/Help_Gene_Ontology>`__ for more details).

example: 

::

	thaliana:
    	star_index: "--genomeSAindexNbases 12"
    	genomesize: 1.3e8
    	ncbiID: "3702" 		# optional
    	genus: "Arabidopsis"	# optional
    	go_database: org.Athaliana.eg.db 	# optional
  
- Resources allocation
   
Output options
--------------

Default parameters
^^^^^^^^^^^^^^^^^^

- Full analysis:
	By default, a full analysis is performed form raw data to analysis plots. Change ``full_analysis`` in the config file.

- Limited QC output:
	By default, some QC options are not performed to limit the time and amount of output files. Change ``QC_option`` in the config file.

- No Gene Ontology analysis: 
	Due to the difficulty in automating building a GO database, this option is OFF by default. Change ``GO`` option in the config file. Please refer to Additional output options #2 below and `Help GO <https://github.com/joncahn/epigeneticbutton/blob/main/Help/Help_Gene_Ontology>`__ before setting it to ``true`` as it requires 2 other files. These files are available for Arabidopsis thaliana (Tair10 / ColCEN assembly) and Maize B73 (v5 or NAM assembly) in the ``data`` folder.

- No TE analysis: 
	By default, no analysis on transposable elements is performed. Change ``te_analysis`` in the config file.

- For ChIP-seq: 
	The default mapping parameters are bowtie2 ``--end-to-end`` default parameters. Other options are available in the config file ``chip_mapping_option``.

- For sRNA-seq: 
	The default is not based on Netflex v3 library preparation. If your data was made with this kit, an additional deduplication and read trimming is required. To turn it ON, change the ``netflex_v3_deduplication`` in the config file. See `Known issues #3` if you have mixed libraries.

    The default is not to filter structural RNAs prior to shortstack analysis. Change ``structural_rna_depletion`` in the config file.  While this step is recommended for small interfering RNA analysis, it requires a pre-build database of fasta files. Please refer to the `Help Rfam <https://github.com/joncahn/epigeneticbutton/blob/main/Help/Help_structural_RNAs_database_with_Rfam>`__  before setting it to ``true``. This file is available for Maize in the ``data`` folder. 

    The default is to only perform *de novo* micro RNA identification (``--dn_mirna`` argument in ShortStack). If you also want the known microRNAs, download the fasta file from `miRbase <https://www.mirbase.org>`__, filter it for your species of interest, and add to the ``srna_mapping_params`` entry in the config file ``--known_miRNAs <path/to/known_miRNA_file.fa>``.

Configuration Options
^^^^^^^^^^^^^^^^^^^^^

- Main output options

	+ ``full_analysis``: When ``false``, only the mapping and the bigwigs will occur. When ``true``, will also be performed: single-data analyses (e.g. peak calling for ChIP, differential expression for RNAseq, DMRs for mC) and combined analyses (e.g. Upset plots for ChIP/TF, heatmaps and metaplots on all genes).

	+ ``te_analysis``: When ``true``, small RNA differential expression will be performed (if such data is available), as well as heatmaps and metaplots of all the samples. The name and path to the TE file in bed format must be filled in the config file for the corresponding reference genome. The name of the TEs (4th column of the bed file) must be unique.

	+ ``QC_option``: When ``true``, runs fastQC on raw and trimmed fastq files.

- ChIP Mapping Parameters

	+ ``default``: Standard mapping parameters
	+ ``repeat``: Centromere-specific mapping (more sensitive)
	+ ``repeatall``: Centromere mapping with relaxed MAPQ
	+ ``all``: Relaxed mapping parameters

- DMRs parameters

	+ By default, DNA methylation data will be analyzed in all sequence contexts (CG, CHG and CHH, where H = A, T or C). The option for CG-only is under development.
	+ DMRs are called with the R package `DMRcaller <https://www.bioconductor.org/packages/release/bioc/html/DMRcaller.html>`__ (DOI: 10.18129/B9.bioc.DMRcaller) for CG, CHG and CHH and the following (stringent) parameters:

::

		-CG: method="noise-filter", binSize=200, test="score", pValueThreshold=0.01, minCytosinesCount=5, minProportionDifference=0.3, minGap=200, minSize=50, minReadsPerCytosine=3	
		-CHG: method="noise_filter", binSize=200, test="score", pValueThreshold=0.01, minCytosinesCount=5, minProportionDifference=0.2, minGap=200, minSize=50, minReadsPerCytosine=3
		-CHH: method="bins", binSize=200, test="score", pValueThreshold=0.01, minCytosinesCount=5, minProportionDifference=0.1, minGap=200, minSize=50, minReadsPerCytosine=3
	
	+ These parameters were selected based on the most optimal results obtained by the authors `[Catoni et al. 2018] <https://academic.oup.com/nar/article/46/19/e114/5050634>`__.
	+ A deeper analysis is available to try different parameters and methods to call the DMRs. Toggle the ``custom_script_dmrs`` on the ``config/config.yaml`` file to use it, and feel free to edit it as well to try different parameters.

- In-line customization of the parameters

More options can be customized by editing the config file directly, and are relatively self-explanatory. Another option is to add the option to the snakemake command directly with the ``--config`` parameter.

For example:
``snakemake --cores 1 --config chip_mapping_option="repeatall"`` will override the ChIP-seq mapping option in the config file and instead use the "repeatall" option for ChIP-seq data in the run. Any argument given in line will take precedent over the value in the config file.

Resources and Profiles
======================

Summary
-------

Resources are pre-defined based on efficient yet conservative values for optimal time/resources requirements ratio. 

It is recommended to check that the default values are adapted to your system is recommended. It includes two aspects: 
   - the rule-specific resources requirements (at the bottom of the ``config/config.yaml`` file)
   - the profile depending on the cluster manager (e.g. for slurm, in ``profiles/slurm/config.yaml``)

At the moment, the values are adapated to the CSHL cluster.

Rule-specific resources requirements
------------------------------------

Each rule has been assigned a set of system resources based on the computational requirements of the task performed. Feel free to edit each set or the specific rule in the config file if you want to allow more or less resources to be used. Each set has a number of threads to use in parallel, a maximum memory and tmp values (in mb) and a "quality of service" (qos) only used for slurm profiles.
The number of threads used in total (potentially by multiple rules in parallel) is capped by the number defined by the ``--cores`` parameters or in the profile. The 

Default sets of allocated resources:

::

	low_resources:
	  threads: 1
	  mem_mb: 1000
	  tmp_mb: 1000
	  qos: "default"

	standard_resources:
	  threads: 4
	  mem_mb: 2000
	  tmp_mb: 2000
	  qos: "default"

	heavy_resources:
	  threads: 8
	  mem_mb: 16000
	  tmp_mb: 48000
	  qos: "slow_nice"
  
	max_resources:
	  threads: 16
	  mem_mb: 32000
	  tmp_mb: 96000
	  qos: "slow_nice"

	single_thread:
	  threads: 1
	  mem_mb: 32000
	  tmp_mb: 48000
	  qos: "slow_nice"

Profiles
--------

If you are running the pipeline on a different platform than CSHL slurm cluster, you will likely need to adjust the config file for your cluster scheduler (`profiles/slurm/config.yaml` for SLURM or create a new profile for your scheduler). It cannot be easily automatized since each scheduler has often their own vocabulary. Ask your cluster IT support for help to set the correct parameters. You might need to install a specific snakemake-executor-plugin based on the type of system you want to deploy the epigenetic button on. Feel free to ask for help or open an issue if needed. 
