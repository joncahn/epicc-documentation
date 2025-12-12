.. _installation:

Installation
------------

1. Clone the repository:

.. code-block:: console

  git clone https://github.com/joncahn/epigeneticbutton.git

or for ssh connection

.. code-block:: console
  
  git clone git@github.com:joncahn/epigeneticbutton.git
  cd epigeneticbutton

2. Install snakemake and other required packages from the depency file:

.. code-block:: console
  
  conda create -n smk9 -y --file config/smk9.txt
  conda activate smk9

If you don't have conda yet, follow the directions here:
https://docs.conda.io/projects/conda/en/stable/user-guide/install/

If you want to run the pipeline on a different platform than locally or slurm, you will need to also install the corresponding snakemake-executor-plugin
