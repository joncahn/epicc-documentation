Installation
------------

1. Clone the repository:

``git clone https://github.com/joncahn/epigeneticbutton.git``

or for ssh connection

``git clone git@github.com:joncahn/epigeneticbutton.git``

``cd epigeneticbutton``

2. Install snakemake and other required packages from the depency file:

``conda create -n smk9 -y --file config/smk9.txt``

``conda activate smk9``

If you don't have conda yet, follow the directions `on the conda website linked here`_

.. _on the conda website linked here: https://docs.conda.io/projects/conda/en/stable/user-guide/install/

If you want to run the pipeline on a different platform than locally or slurm, you will need to also install the corresponding snakemake-executor-plugin (`executor plugin options here`_)

.. _executor plugin options here: https://snakemake.github.io/snakemake-plugin-catalog/index.html
