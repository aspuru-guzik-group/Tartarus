# Tartarus: Practical and Realistic Benchmarks for Inverse Molecular Design

This repository contains the code and results for the paper [Tartarus, an open-source collection of benchmarks for evaluation of a generative model](https://arxiv.org/abs/2209.12487).

Total installation time: ~15-20mins. 

## Installing XTB and CREST

The task of designing organic photovoltaics and emitters will require the use of [**XTB**](https://github.com/grimme-lab/xtb), a program package of semi-empirical quantum mechanical methods, and [**CREST**](https://github.com/grimme-lab/crest), a utility of xtb used to sample molecular conformers. 

The binaries are provided in [here](https://github.com/grimme-lab/xtb/releases). Place in home directory, and the software can be sourced using
```bash
export XTBHOME=${HOME}/xtb
export PATH=${PATH}:${XTBHOME}/bin
export XTBPATH=${XTBHOME}/share/xtb:${XTBHOME}:${HOME}
export MANPATH=${MANPATH}:${XTBHOME}/share/man
```

## Installing SMINA

The task of designing molecules that dock to proteins requires the use of [**SMINA**](https://sourceforge.net/projects/smina/), a method for calcualte docking scores of ligands onto solved structures (proteins). The binary file is already included in the repository, in `tartarus/docking_structures/smina.static`.


## Packages required
Use `python >= 3.8`. We recommend using a conda environment for the installation of 
- rdkit >= 2021.03.3
- xtb-python >= 20.1
- openbabel == 3.1.1
<!-- - xtb == 6.4.1
- xorg-libxrender == 0.9.10 -->

Required packages:
- numpy >= 1.22.3
- pandas >= 1.4.3 
- torch == 1.12.0
- pyscf == 2.0.1
- morfeus-ml >= 0.7.1
- geometric == 0.9.7.2
- pyberny == 0.6.3
- loguru == 0.6.0
- geodesic-interpolate == 1.0.0 \
(`pip install -i https://test.pypi.org/simple/ geodesic-interpolate`)
- polanyi == 0.0.1 \
(`pip install git+https://github.com/kjelljorner/polanyi`)
<!-- - selfies == 1.0.3
- scikit-learn >= 1.1.1 
- h5py == 3.7.0
- wurlitzer == 3.0.2
- sqlalchemy >= 1.4.13 -->


## Datasets 

All datasets are found in the [datasets](datasets/) directory. The arrows indicate the goal (&#8593; = maximization, &#8595; = minimization). 

|Task | Dataset name       | # of smiles |  Columns in file ||||||
|---|--------------------|------------------|----|----|----|---|----|----|
| Designing OPV | `hce.csv`          | 24,953         | PCE<sub>PCBM</sub> -SAS (&#8593;) | PCE<sub>PCDTBT</sub> -SAS (&#8593;) | 
| Designing emitters | `gdb13.csv`        | 403,947          | Singlet-triplet gap (&#8595;) | Oscillator strength (&#8593;) | Multi-objective (&#8593;) |  | ||
| Designing drugs | `docking.csv`      | 152,296          | 1SYH (&#8595;) | 6Y2F (&#8595;) | 4LDE (&#8595;) |  | | |
| Designing chemical reaction substrates | `reactivity.csv`      | 60,828          | 	Activation energy &Delta;E<sup>&#8225;</sup> (&#8595;)   |  	Reaction energy &Delta;E<sub>r</sub> (&#8595;)  | &Delta;E<sup>&#8225;</sup> + &Delta;E<sub>r</sub> (&#8595;)  |  - &Delta;E<sup>&#8225;</sup> + &Delta;E<sub>r</sub> (&#8595;)    |     | |

## Getting started 

Below are some examples of how to load the datasets and use the fitness functions. For more details, you can also look at `example.py`. 

### Designing organic photovoltaics

To use the evaluation function, load either the full xtb calculation from the `pce` module, or use the surrogate model, with pretrained weights.

```python
import pandas as pd
data = pd.read_csv('./datasets/hce.csv')   # or ./dataset/unbiased_hce.csv
smiles = data['smiles'].tolist()
smi = smiles[0]

## use full xtb calculation in hce module
from tartarus import pce
dipm, gap, lumo, combined, pce_pcbm_sas, pce_pcdtbt_sas = pce.get_properties(smi)

## use pretrained surrogate model
dipm, gap, lumo, combined = pce.get_surrogate_properties(smi)
```


### Designing Organic Emitters
Load the objective functions from the `tadf` module. All 3 fitness functions are returned for each smiles.

```python
import pandas as pd
data = pd.read_csv('./datasets/gdb13.csv')  
smiles = data['smiles'].tolist()
smi = smiles[0]

## use full xtb calculation in hce module
from tartarus import tadf
st, osc, combined = tadf.get_properties(smi)
```


### Design of drug molecule
Load the `docking` module. There are separate functions for each of the proteins, as shown below.

```python
import pandas as pd
data = pd.read_csv('./datasets/docking.csv')  
smiles = data['smiles'].tolist()
smi = smiles[0]

## Design of Protein Ligands 
from tartarus import docking
score_1syh = docking.get_1syh_score(smi)
score_6y2f = docking.get_6y2f_score(smi)
score_4lde = docking.get_4lde_score(smi)
```


### Design of Chemical Reaction Substrates
Load the `reactivity` module. All 4 fitness functions are returned for each smiles.

```python
import pandas as pd
data = pd.read_csv('./datasets/reactivity.csv')  
smiles = data['smiles'].tolist()
smi = smiles[0]

## calculating binding affinity for each protein
from tartarus import reactivity
Ea, Er, sum_Ea_Er, diff_Ea_Er = reactivity.get_properties(smi)
```


### Results
Our results for running the corresponding benchmarks can be found here: 
- Design of Protein Ligands: https://drive.google.com/file/d/1d_4mg1Eb7HrUJ2L7A8kFtld-TmPmOKlJ/view?usp=sharing
- Design of Chemical Reaction Substrates: https://drive.google.com/file/d/1fCnFxSUITg4qSlOuwFolvQPUQA31Qaii/view?usp=sharing
- Designing organic photovoltaics (photovoltaic conversion efficiency): https://drive.google.com/file/d/1w6oOBGjDC4Enh492jLQ7A3Xc1XbHXiIt/view?usp=sharing
- Designing Organic Emitters: https://drive.google.com/file/d/1l8weYg835HDGvOoRbOcHUnvLjiyQi_Ms/view?usp=sharing
- Designing organic photovoltaics (Explore): https://drive.google.com/file/d/1-J99iXfBx0_aG1BqEEXPh7q0kovBFD0L/view?usp=sharing
- Designing organic photovoltaics (Surrogate, exploit): https://drive.google.com/file/d/1EV7ST9_F4DBnQpxhd6VaaJWP5r9ygr0c/view?usp=sharing
- Designing organic photovoltaics (Exploit): https://drive.google.com/file/d/1Yh_8E3jRf6X230CvlRlPtk2qPQIkC5hB/view?usp=sharing


## Questions, problems?
Make a github issue 😄. Please be as clear and descriptive as possible. Please feel free to reach
out in person: (akshat98[AT]stanford[DOT]edu, robert[DOT]pollice[AT]gmail[DOT]com)

## License

[Apache License 2.0](https://choosealicense.com/licenses/apache-2.0/)
