====
Help
====

Specific help for additional analysis.

  
Help Structural RNA database with Rfam
======================================

How to get the file of structural RNAs from Rfam database to filter from sRNA-seq in the EPICC pipeline.
Help: `<https://docs.rfam.org/en/latest/sequence-extraction.html>`

- STEP 1: install easel::

    conda install easel

- STEP 2: download all Rfam fasta files from the database::

    wget ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/fasta_files/Rfam.fa.gz .

- STEP 3: unzip and combine into one file::

    gunzip *.gz

- STEP 4: index the fata file::

    esl-sfetch --index Rfam.fa

If the indexing does not work because of repetitive fasta sequences, try::

  awk 'BEGIN{RS=">"; FS="\n"; ORS=""} (FNR==1){next} { name=$1; seq=$0; gsub(/(^[^\n]*|)\n/,"",seq) } !(seen[name,seq]++){ print ">" $0 }' Rfam.fa

- STEP 5: create an `.sql` file (for example `query_rRNA.sql`) that fetches the regions of interest which contains the following code::

    SELECT concat(fr.rfamseq_acc,'/',seq_start,'-',seq_end)
    FROM full_region fr, rfamseq rf, taxonomy tx, family f
    WHERE
    rf.ncbi_id = tx.ncbi_id
    AND f.rfam_acc = fr.rfam_acc
    AND fr.rfamseq_acc = rf.rfamseq_acc
    AND tx.ncbi_id = 4577 ## ncbi ID of the organism or interest, this is Zea mays
    AND f.type LIKE '%rRNA%' ## type of features to fetch
    AND is_significant = 1;

- STEP 6: repeat STEP5 to create .sql files for all type of structural RNAs you want to filter, changing the `f.type` pattern (`%rRNA%` in the example above).

- STEP 7: Get the corresponding accessions from mysql, replacing query.sql by the name of the previously created files::

    mysql -urfamro -hmysql-rfam-public.ebi.ac.uk -P4497 --skip-column-names --database Rfam < query.sql > accessions.txt

If using several files, run this command each time and concatenate the outputs, making sure each entry is unique, and all the accessions are present in the Rfam.fa::

  cat accessions_*.txt | sort -u > unique_accessions.txt
  while read accession
  do
          if [ $(grep $accession Rfam.fa | wc -l) -gt 0 ]; then
                  printf "$accession\n" >> good_accessions.txt
          fi
  done < unique_accessions.txt

- STEP 8: Extract the fasta sequences of the chosen accessions from the Rfam database, replacing the `/path/to/Rfam.fa` and `/path/to/accessions` by the actual path to the downloaded file in STEP1 and the `good_accessions.txt` file from STEP7::

    esl-sfetch -f path/to/Rfam.fa /path/to/accessions.txt > Rfam_ncRNAs.fa

- STEP 9: gzip the file::
  
    gzip Rfam_ncRNAs.fa

- STEP 10: Edit the `structural_rna_fafile` entry in the config file with this file (including the full path).

Gene Ontology analysis
======================

How to create the GO database for use in the EPICC pipeline.

Prexisting database
+++++++++++++++++++

- If you already have a database from AnnotationForge to use, copy the file in the `genomes/<ref_genome>/GO/` directory for the analysis to work automatically. If you have not run the epigenetic button yet, you will need to create these folders, where `ref_genome` is the last column of your sample file.

- For some species, an annotation package can be created with AnnotationForge makeOrgPackageFromNCBI(). For help:
https://bioconductor.org/packages/release/bioc/vignettes/AnnotationForge/inst/doc/MakingNewOrganismPackages.html

Creating a new database
+++++++++++++++++++++++

To create a new database *de novo*, two files are required, a file linking genes to GO (``gaf_file``: <ref>_infoGO.tab, originally GAF format) and a file with information about the genes (``gene_info_file``: <ref>_genes_info.tab).

These steps are to be performed before the analysis run for the GO analysis to be performed automatically with the rest of the analysis.
If created after, run snakemake with the GO file of interest as a target. See `:ref:<>`

- STEP 1: Get the GAF file (``gaf_file``)

The GAF file can usually be downloaded from NCBI, or from species-specific community resources (e.g. https://www.arabidopsis.org/ for Arabidopsis).

Help to find the GAF file: https://geneontology.org/docs/download-go-annotations/ 

Importantly, this file should be tab separated have gene IDs in the first column (e.g. AT1G00010), the GO terms in the 6th column (e.g GO:00001) and the evidence in the 10th column (i.e. IEA). The other columns do not matter and can be empty.

If you do not want to edit the file itself, you can instead edit the scripts ``R_build_GO_database.R`` and ``R_GO_analysis`` in ``workflow/scripts/`` for the fGO table to select the right columns.

Add this file to the corresponding reference genome ``gaf_file`` in the config file::

  B73_v5:
    gaf_file: data/B73_v5_infoGO.tab

- STEP 2: Create the gene information file (``gene_info_file``)

The ``gene_info_file`` can be easily generated from a gff file.

The fourth column ``GID`` must match the ``GID`` from the first column of the ``gaf_file`` above. 
 
For example (for B73_v5)::

  printf "Chr\tStart\tEnd\tGID\tType\tDescription\n" > data/B73_v5_genes_info.tab
  awk -v OFS="\t" '$3=="gene" {print $1,$4-1,$5,$9,".",$7}' genomes/B73_v5/B73_v5.gff | awk -F"[:;=]" -v OFS="\t" '{print $1,$2,$4,$6}' | awk -v OFS="\t" '{print   $1,$2,$3,$5,$6,$7}' >> data/B73_v5_genes_info.tab

Add this file to the corresponding reference genome ``gene_info_file`` in the config file::

  B73_v5:
    gene_info_file: data/B73_v5_genes_info.tab

- STEP 4: Run the epigenetic button once only to generate the GO database::

    snakemake --cores 1 genomes/<ref_genome>/GO/<dbname> 

for example::

  snakemake --cores 1 genomes/ColCEN/GO/org.Zmays.eg.db

- STEP 4b: Create the GO database manually

Run the script: ``workflow/scripts/R_build_GO_database.R``, with the following arguments::

  script="scripts/R_build_GO_database.R"
  infofile="B73_v5_infoGO.tab" # replace with appopriate `gaf_file` from STEP1
  genefile="B73_v5_genes_info.tab # replace with appropriate `gene_info_file` from STEP2
  ref_genome="B73_v5" # replace with corresponding genome reference (same than on the sample file)
  genus="Zea" # replace with corresponding genus (conventionally first letter is capitalized)
  species="mays" # replace with corresponding species (conventionally lowerscript)
  ncbiID="4577" # find value corresponding to the species at NCBI

  Rscript ${script} ${infofile} ${genefile} ${ref_genome} ${genus} ${species} ${ncbiID}

An replace the `go_database` entry in the config file for the corresponding species:: 

  mays:
    go_database: org.Zmays.eg.db # where dbanme="org.<firstlettergenus><species>.eg.db"
                                                                                                                        
                                                                                                                        
                                                                                                                        
