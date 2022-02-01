import rdkit
from rdkit.Chem import AllChem
from rdkit.Chem import AllChem as Chem
from torch.utils.data import Dataset, DataLoader
from rdkit.Chem import DataStructs
import numpy as np
import os, sys
#sys.path.append(os.path.abspath(os.path.join(__file__, "..", "..", "..", "..")))
#from autoPyTorch import AutoNetRegression
import pandas as pd
import json
import matplotlib.pyplot as plt
import time
import torch
import torch.nn as nn
from sklearn.metrics import r2_score

def get_fingerprint(smile, nBits, ecfp_degree=2):
    m1 = Chem.MolFromSmiles(smile)
    fp = AllChem.GetMorganFingerprintAsBitVect(m1,ecfp_degree, nBits=nBits)
    x = np.zeros((0, ), dtype=np.int8)
    DataStructs.ConvertToNumpyArray(fp, x)

    return x


def load_hce_data(suffix):
    X = []
    Y_homo_lumo = []
    Y_lumo_val = []
    Y_dipole = []
    Y_function = []

    with open('./hce_smiles_{}.smi'.format(suffix), 'r') as smiles:
        for smi in smiles.readlines():
            bit_vector = get_fingerprint(smile=smi, nBits=2048, ecfp_degree=2)
            X.append(bit_vector)
        
    with open('./hce_smiles_props_{}.txt'.format(suffix), 'r') as props:
        for prop in props.readlines():
            prop = prop.split(' ')
            homo_lumo, lumo_val, dipole, function = prop[0], prop[1], prop[2], prop[3]
            Y_homo_lumo.append(float(homo_lumo))
            Y_lumo_val.append(float(lumo_val))
            Y_dipole.append(float(dipole))
            Y_function.append(float(function))

    return np.array(X, dtype=np.float32), np.array(Y_homo_lumo, dtype=np.float32), np.array(Y_lumo_val, dtype=np.float32), np.array(Y_dipole, dtype=np.float32), np.array(Y_function, dtype=np.float32)



class SurrogateModel(object):
    def __init__(self, model_list_dir, use_ensemble=False):
        super(SurrogateModel, self).__init__()
        self.use_ensemble = use_ensemble
        
        model_state_dicts = os.listdir(model_list_dir)
        self.model_list = []
        for model_state_dict in model_state_dicts:
            self.model_list.append(torch.load(os.path.join(model_list_dir, model_state_dict)))

            print(torch.load(os.path.join(model_list_dir, model_state_dict)))

            if use_ensemble is False:
                break

        

    def forward(self, x):
        if type(x) == str:
            x = get_fingerprint(smile=x, nBits=2048, ecfp_degree=2)
            
        x = torch.tensor(x).to('cuda:0', dtype=torch.float32)
        predictions = []
        for model in self.model_list:
            predictions.append(model(x).detach().cpu().numpy())

        predictions = np.array(predictions)

        mean = np.mean(predictions, axis=0)
        var = np.var(predictions, axis=0)

        if self.use_ensemble:
            return mean, var
        else:
            return mean

        
def calc_properties_with_surrogate(smiles):
    try: 
        surrogate_model_homo_lumo = SurrogateModel('./trained_model/homo_lumo')
        surrogate_model_lumo_val = SurrogateModel('./trained_model/lumo_val')
        surrogate_model_dipole = SurrogateModel('./trained_model/dipole')
        surrogate_model_function = SurrogateModel('./trained_model/function')
        
        homo_lumo = surrogate_model_homo_lumo.forward(smiles) 
        lumo_val = surrogate_model_lumo_val.forward(smiles)
        dipole = surrogate_model_dipole.forward(smiles)
        function = surrogate_model_function.forward(smiles)
    
        return float(dipole), float(homo_lumo), float(-lumo_val), float(function)
    
    except: 
        return float(-10**4), float(-10**4), float(-10**4), float(-10**4)




if __name__ == '__main__':
    
    smi = '[SiH2]1C=Cc2c1c1c3c[nH]cc3c3cc([se]c3c1c1cscc21)-c1scc2sccc12'
    A = calc_properties_with_surrogate(smi)


