Installation
------------

1. Clone the repository:

.. code-block:: bash

   git clone https://github.com/joncahn/epigeneticbutton.git

To clone and rename the directory (recommended):

.. code-block:: bash

   git clone https://github.com/joncahn/epigeneticbutton.git epicc
   cd epicc

Or via SSH:

.. code-block:: bash

   git clone git@github.com:joncahn/epigeneticbutton.git epicc
   cd epicc

2. Install Snakemake and other required packages from the dependency file and activate the environment:

.. code-block:: bash

   conda create -n epicc -y --file config/epicc-env.txt
   conda activate epicc

If you do not have conda yet, follow the directions `on the conda website linked here`_.

.. _on the conda website linked here: https://docs.conda.io/projects/conda/en/stable/user-guide/install/

If you want to run the pipeline on a platform other than a local workstation or a SLURM cluster, you will also need to install the corresponding Snakemake executor plugin (`executor plugin options here`_).

.. _executor plugin options here: https://snakemake.github.io/snakemake-plugin-catalog/index.html
