Configuration
=============

For new users, it is recommended to use the configuration app to validate your sample metadata file and choose analysis options:
https://epicc-builder.streamlit.app/

1. Prepare your sample metadata file (default to `config/all_samples.tsv`) with the required columns below (see Input requirements for more details specific to each data-type):
   
  - `data_type`: Type of data [RNAseq | ChIP_* | TF_* | mC | sRNA] (RAMPAGE under development)
  - `line`: Sample line (e.g. B73)
  - `tissue`: Tissue type
  - `sample_type`: Sample identifier
  - `replicate`: Replicate ID
  - `seq_id`: Sequencing ID; use the corresponding SRR####### if downloading from SRA
  - `fastq_path`: Path to FASTQ files; if downloading from SRA, use "SRA" 
  - `paired`: [PE | SE]
  - `ref_genome`: Reference genome name

2. Update `config/config.yaml` with your paths and parameters:
  
  - Sample file: this is the full path to the file detailed above which contain your samples metadata. 
  - Reference genome files: for each reference genome in the sample file (last column), enter the full path of a fasta file, a gene gff file, and a gene gtf file (See [below](#common-to-all-types-of-samples) for more details)
  - Analysis parameters / options
  - Species-specific parameters
  - Resources allocation
   
3. If you are running the pipeline on a different platform than CSHL slurm cluster, you will likely need to adjust the rule-specific resource parameters at the bottom of the `config/config.yaml` and the config file for your cluster scheduler (`profiles/slurm/config.yaml` for SLURM or create a new profile for your scheduler). In slurm, the default is to start 16 jobs maximum in parallel. Keep in mind that units in the cluster file are in MB.

4. Default options:

  - Full analysis: By default, a full analysis is performed form raw data to analysis plots. Change `full_analysis` in the config file ([see below](#main-output-options)).
  - Limited QC output: By default, some QC options are not performed to limit the time and amount of output files. Change `QC_option` in the config file ([see below](#main-output-options)).
  - No Gene Ontology analysis: Due to the difficulty in automating building a GO database, this option is OFF by default. Change `GO` option in the config file. Please refer to Additional output options #2 below and [Help GO](Help/Help_Gene_Ontology) before setting it to `true` as it requires 2 other files. These files are available for Arabidopsis thaliana (Tair10 / ColCEN assembly) and Maize B73 (v5 or NAM assembly) in the `data` folder.
  - No TE analysis: By default, no analysis on transposable elements is performed. Change `te_analysis` in the config file ([see below](#main-output-options)).
  - For ChIP-seq: the default mapping parameters are bowtie2 `--end-to-end` default parameters. Other options are available in the config file `chip_mapping_option` ([see below](#chip-mapping-parameters)).
  - For sRNA-seq: the default is not based on Netflex v3 library preparation. If your data was made with this kit, an additional deduplication and read trimming is required. To turn it ON, change the `Netflex_v3_deduplication` in the config file. See [Known issues #3](#known-potential-issues) below if you have mixed libraries.
  - For sRNA-seq: the default is not to filter structural RNAs prior to shortstack analysis. Change `structural_rna_depletion` in the config file.  While this step is recommended for small interfering RNA analysis, it requires a pre-build database of fasta files. Please refer to the [Help structural RNAs](Help/Help_structural_RNAs_database_with_Rfam) before setting it to `true`. This file is available for Maize in the `data` folder.
  - For sRNA-seq: the default is to only perform *de novo* micro RNA identification (`--dn_mirna` argument in ShortStack). If you also want the known microRNAs, download the fasta file from [miRbase](https://www.mirbase.org), filter it for your species of interest, and add to the `srna_mapping_params` entry in the config file `--known_miRNAs <path/to/known_miRNA_file.fa>`.
