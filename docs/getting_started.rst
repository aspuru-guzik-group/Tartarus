Getting Started
===============

Install and benchmark with Tartarus quickly!

Setup & Installation
----------------------

There are two options to setup Tartarus. The first is to use the provided Docker image, the second is to install Tartarus on your local machine.

Docker Setup
************

To use the Docker image, you will need to have Docker installed on your machine. You can find instructions for installing Docker `here <https://docs.docker.com/get-docker/>`_.

1. Once Docker is installed, you can pull the Tartarus image from Docker Hub:

.. code-block:: console

    $ docker pull aspuru/tartarus:latest

2. You can then run the image with:

.. code-block:: console

    $ docker run -it aspuru/tartarus:latest --task <task_name> --input_filename <input_filename>

Local Installation
******************

To install Tartarus locally, you will need to have Conda installed on your machine. You can find instructions for installing Conda `here <https://docs.conda.io/projects/conda/en/latest/user-guide/install/>`_. 

1. Clone the Tartarus repository.

.. code-block:: console

    $ git clone git@github.com:aspuru-guzik-group/Tartarus.git

2. Create a Conda environment using the provided environment configuration.

.. code-block:: console

    $ conda env create -f environment.yml

3. Activate the tartarus Conda environment.

.. code-block:: console

    $ conda activate tartarus

4. Set environment variables.

.. code-block:: console

    $ export XTBHOME=$CONDA_PREFIX
    $ source $CONDA_PREFIX/share/xtb/config_env.bash

5. Optionally, you can can configure the environment variables to be set automatically when you activate the Conda environment.

.. code-block:: console

    $ echo "export XTBHOME=$CONDA_PREFIX" > $CONDA_PREFIX/etc/conda/activate.d/env.sh
    $ echo "source $CONDA_PREFIX/share/xtb/config_env.bash" >> $CONDA_PREFIX/etc/conda/activate.d/env.sh

Common Issues
*************


