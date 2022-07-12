#!/bin/bash
conda create --name tartarus python=3.8 -y
conda activate tartarus

# should be installed already
conda install -c conda-forge rdkit openbabel -y
conda install -c pytorch pytorch -y     # install using conda to avoid conflict with libxtb

conda install -c conda-forge xtb-python -y
conda install -c conda-forge crest -y
export XTBHOME=$CONDA_PREFIX
source $CONDA_PREFIX/share/xtb/config_env.bash

pip install --upgrade pip
pip install numpy
pip install pyscf morfeus-ml

# additional packages for polanyi
pip install h5py scikit-learn geometric pyberny loguru wurlitzer sqlalchemy
pip install -i https://test.pypi.org/simple/ geodesic-interpolate
pip install git+https://github.com/kjelljorner/polanyi

