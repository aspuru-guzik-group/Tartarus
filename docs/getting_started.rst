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

    $ docker pull johnwilles/tartarus:latest

2. You can then run the image with:

.. code-block:: console

    $ docker run -v /local/path/to/data:/data johnwilles/tartarus:latest --mode <mode_name> --input_filename <input_filename>

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

6. Ensure that docking task executables have the correct permissions.

.. code-block:: console

    $ chmod 777 tartarus/data/qvina
    $ chmod 777 tartarus/data/smina

Common Issues
*************

Depending on the version of Conda that you have installed, it is possible that the geodesic-interpolate package may not install correctly 
from the PyPI test registery. If this is the case, you can install the package manually by running:

.. code-block:: console

    $ pip install --extra-index-url https://test.pypi.org/simple/ geodesic-interpolate

Benchmarking Quick Start
------------------------

The quickest way to get started with Tartarus is to use the provided Docker image. You can run the image with the following command:

.. code-block:: console

    $ docker run -v /local/path/to/data:/data johnwilles/tartarus:latest --mode <mode_name> --input_filename <input_filename>

The Docker ``-v`` flag mounts the local directory ``/local/path/to/data`` to the Docker container's ``/data`` directory. This allows Tartarus to access the data inside the container. The benchmarking script exposes the following configuration flags:

* ``--mode``: The name of the benchmarking mode to run. The available modes are ``pce``, ``tadf``, ``docking``, ``reactivity``.
* ``--input_filename``: The name of the input file to use for the benchmarking task. The input file should be located in the mounted directory.
* ``--output_filename``: The name of the output file to write the benchmarking results to. The output file will be written to the mounted directory. If this flag is not provided, the results will be written to ``output.csv``.
* ``--parallel``: Configures Tartarus to use parallel processes for the benchmarking task. If this flag is provided, Tartarus will use all available cores. If this flag is not provided, Tartarus will use a single core.
* ``--verbose``: Configures Tartarus to print verbose output to the console. If this flag is provided, Tartarus will print verbose output.