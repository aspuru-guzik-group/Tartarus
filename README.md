# Tartarus: The Next Generation of Benchmarks for Inverse Molecular Design

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

The task of designing molecules that dock to proteins requires the use of [**SMINA**](https://sourceforge.net/projects/smina/), a method for calcualte docking scores of ligands onto solved structures (proteins).


## Datasets 

All datasets are found in the [datasets](datasets/) directory. 

|Task | Dataset name       | Number of smiles |  Columns in file |||||
|---|--------------------|------------------|----|----|----|--|---|
| Designing OPV | `hce.csv`          | 24,953           | HOMO-LUMO gap (&#8593;) | LUMO (&#8595;) | Dipole (&#8593;) | Combined objective (&#8593;) |
| Designing OPV | `unbiased_hce.csv` | 1,000            | HOMO-LUMO gap (&#8593;) | LUMO (&#8595;) | Dipole (&#8593;) | Combined objective (&#8593;) |
| Designing emitters | `gdb13.csv`        | 403,947          | Singlet-triplet gap (&#8595;) | Oscillator strength (&#8593;) | Multi-objective (&#8593;) | Time (s) |
| Designing drugs | `docking.csv`      | 152,296          | 1SYH (&#8595;) | 6Y2F (&#8593;) | 4LDE (&#8593;) | Time (s) | |
| Designing drugs | `reactivity.csv`      | 60,828          | |  | |  | |


## Getting started 


### Designing organic photovoltaics

To use the evaluation function, load either the full xtb calculation from the `pce` module, or use the surrogate model, with pretrained weights.

```python
import pandas as pd
data = pd.read_csv('./datasets/hce.csv')   # or ./dataset/unbiased_hce.csv
smiles = data['smiles'].tolist()
smi = smiles[0]

## use full xtb calculation in hce module
from tartarus import pce
dipole, hl_gap, lumo, combined, pce_1, pce_2, sas = pce.get_properties(smi)

## use pretrained surrogate model
dipole, hl_gap, lumo, combined = pce.get_surrogate_properties(smi)
```


### Designing Organic Emitters

```python
import pandas as pd
data = pd.read_csv('./datasets/gdb13.csv')   # or ./dataset/unbiased_hce.csv
smiles = data['smiles'].tolist()
smi = smiles[0]

## use full xtb calculation in hce module
from tartarus import tadf
st, osc, combined = tadf.get_properties(smi)
```


### Design of drug molecule

```python
import pandas as pd
data = pd.read_csv('./datasets/docking.csv')   # or ./dataset/unbiased_hce.csv
smiles = data['smiles'].tolist()
smi = smiles[0]

## Design of Protein Ligands 
from tartarus import docking
st, osc, combined = docking.get_1syh_score(smi)
st, osc, combined = docking.get_6y2f_score(smi)
st, osc, combined = docking.get_4lde_score(smi)
```


### Design of Chemical Reaction Substrates

```python
import pandas as pd
data = pd.read_csv('./datasets/reactivity.csv')   # or ./dataset/unbiased_hce.csv
smiles = data['smiles'].tolist()
smi = smiles[0]

## calculating binding affinity for each protein
from tartarus import reactivity
activation_energy, reaction_energy, sa_score = reactivity.get_properties(smi)
```


### Results
Our results for running the corresponding benchmarks can be found here: 
- Design of Protein Ligands: https://drive.google.com/file/d/1d_4mg1Eb7HrUJ2L7A8kFtld-TmPmOKlJ/view?usp=sharing









