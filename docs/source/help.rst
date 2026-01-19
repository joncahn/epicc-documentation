Help
====

Specific help for additional analysis.

  
Help Structural RNA database with Rfam
======================================

Downloading the structural RNAs files from Rfam database to filter from sRNA-seq
Help: https://docs.rfam.org/en/latest/sequence-extraction.html

-STEP1: install easel (conda install easel)

-STEP2: download all Rfam fasta files from the database

wget ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/fasta_files/Rfam.fa.gz .

-STEP3: unzip and combine into one file

gunzip *.gz

-STEP4: index the fata file

esl-sfetch --index Rfam.fa

If the indexing does not work because of repetitive fasta sequences, try:

awk 'BEGIN{RS=">"; FS="\n"; ORS=""} (FNR==1){next} { name=$1; seq=$0; gsub(/(^[^\n]*|)\n/,"",seq) } !(seen[name,seq]++){ print ">" $0 }' Rfam.fa

-STEP5: create the .sql files that fetches the regions of interest

SELECT concat(fr.rfamseq_acc,'/',seq_start,'-',seq_end)
FROM full_region fr, rfamseq rf, taxonomy tx, family f
WHERE
rf.ncbi_id = tx.ncbi_id
AND f.rfam_acc = fr.rfam_acc
AND fr.rfamseq_acc = rf.rfamseq_acc
AND tx.ncbi_id = 4577 ## ncbi ID of the organism or interest, this is Zea mays
AND f.type LIKE '%rRNA%' ## type of features to fetch
AND is_significant = 1;

-STEP6: Create .sql files for all type of structural RNAs you want to filter together

- STEP7: Get the corresponding accessions from mysql, replacing query.sql by the name of the previously created files.

mysql -urfamro -hmysql-rfam-public.ebi.ac.uk -P4497 --skip-column-names --database Rfam < query.sql > accessions.txt
If using several files, run this command each time and concatenate the outputs, making sure each entry is unique, and all the accessions are present in the Rfam.fa:
cat accessions_*.txt | sort -u > unique_accessions.txt
while read accession
do
        if [ $(grep $accession Rfam.fa | wc -l) -gt 0 ]; then
                printf "$accession\n" >> good_accessions.txt
        fi
done < unique_accessions.txt

- STEP8: Extract the fasta sequences of the chosen accessions from the Rfam database

esl-sfetch -f path/to/Rfam.fa /path/to/accessions.txt > Rfam_ncRNAs.fa

- STEP9: gzip and put the path to this file in the config file <structural_rna_fafile>

                                                                                                                        
Gene Ontology analysis
======================
                                                                                                                        
## To create the GO database

If you already have a database from AnnotationForge to use, copy the file in the "genomes/<ref_genome>/GO/" directory for the analysis to work automatically. 
For some species, an annotation package can be created with AnnotationForge makeOrgPackageFromNCBI(). For help:
https://bioconductor.org/packages/release/bioc/vignettes/AnnotationForge/inst/doc/MakingNewOrganismPackages.html

To create a new package de novo, two files are required, a file linking genes to GO (FILE1=<ref>_infoGO.tab, originally GAF format) and a file with information about the genes (FILE2=<ref>_genes_info.tab)
The GAF file can usually be downloaded from NCBI, or from species-specific community resources (e.g. https://www.arabidopsis.org/ for Arabidopsis).
Help to find the GAF file: https://geneontology.org/docs/download-go-annotations/ 
The basic gene information file can be generated from a gff file.

Example for B73 (NAM/v5 version):

File1:
awk '$1 !~ /^!/' maize.B73.AGPv2.aggregate.gaf.gz > data/B73_v5_infoGO.tab 
# Importantly, this file should have gene IDs in the first column (e.g. AT1G00010), the GO terms in the 6th column (e.g GO:00001) and the evidence in the 10th column (i.e. IEA)
# Otherwise, edit the R_build_GO_database.R for the fGO table to select the right columns, as well as the R_GO_analysis.

File2:
printf "Chr\tStart\tEnd\tGID\tType\tDescription\n" > data/B73_v5_genes_info.tab
awk -v OFS="\t" '$3=="gene" {print $1,$4-1,$5,$9,".",$7}' genomes/B73_v5/B73_v5.gff | awk -F"[:;=]" -v OFS="\t" '{print $1,$2,$4,$6}' | awk -v OFS="\t" '{print $1,$2,$3,$5,$6,$7}' >> data/B73_v5_genes_info.tab
The first column GID must match the GID from the first column of FILE1. 

Then, run the script: R_build_GO_database.R, with the following arguments:
script="scripts/R_build_GO_database.R"
ïnfofile="B73_v5_infoGO.tab" # replace with modified FILE1
genefile="B73_v5_genes_info.tab # replace with modified FILE2
ref_genome="B73_v5" # replace with corresponding genome reference (same than on the sample file)
genus="Zea" # replace with corresponding genus (conventionally first letter is capitalized)
species="mays" # replace with corresponding species (conventionally lowerscript)
ncbiID="4577" # find correpsonding value at NCBI

Rscript ${script} ${infofile} ${genefile} ${ref_genome} ${genus} ${species} ${ncbiID}

and replace the GOdatabase in the config file of the corresponding species with: "org.Zmays.eg.db", where dbanme="org.<firstlettergenus><species>.eg.db"

Or fill the config file appropriately and run:
snakemake --cores 1 genomes/<ref_genome>/GO/<dbname> 
for example:
snakemake --cores 1 genomes/ColCEN/GO/org.Zmays.eg.db

These steps are to be performed before the analysis run for the GO analysis to be performed automatically with the rest of the analysis.
If created after, run snakemake with the GO file of interest as a target.
                                                                                                                        
                                                                                                                        
                                                                                                                        
