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
