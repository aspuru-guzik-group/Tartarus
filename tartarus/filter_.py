# -*- coding: utf-8 -*-
"""
Created on Fri Jun  9 11:58:01 2023

@author: aksha
"""
import os 
import pandas as pd 
from rdkit.Chem import MolFromSmiles as smi2mol
from rdkit import Chem
from rdkit import RDLogger
from rdkit.Chem.Crippen import MolLogP
from rdkit.Chem.Descriptors import ExactMolWt
from rdkit.Chem.Lipinski import NumHAcceptors, NumHDonors
import rdkit.Chem as rdc
import rdkit.Chem.rdmolops as rdcmo
import rdkit.Chem.Descriptors as rdcd
import rdkit.Chem.rdMolDescriptors as rdcmd
from rdkit.Chem import QED
from rdkit.Chem import Descriptors
from rdkit.Chem import RDConfig
import sys
sys.path.append(os.path.join(RDConfig.RDContribDir, 'SA_Score'))
import sascorer

import os    
RDLogger.DisableLog('rdApp.*')


_mcf = pd.read_csv(os.path.join('./data/mcf.csv'))
_pains = pd.read_csv(os.path.join('./data/wehi_pains.csv'), names=['smarts', 'names'])
_filters = [Chem.MolFromSmarts(x) for x in _mcf.append(_pains, sort=True)['smarts'].values]


inf = open("./data/pains.txt", "r")
sub_strct = [ line.rstrip().split(" ") for line in inf ]
smarts = [ line[0] for line in sub_strct]
desc = [ line[1] for line in sub_strct]
dic = dict(zip(smarts, desc))


def lipinski_filter(smiles):
    mol = Chem.MolFromSmiles(smiles)
    try: 
        return MolLogP(mol) <= 5 and NumHAcceptors(mol) <= 10 and NumHDonors(mol) <= 5 and 100 <= ExactMolWt(mol) <= 500
    except: 
        return False
    
    
def maximum_ring_size(mol):
    """
    Calculate maximum ring size of molecule
    """
    cycles = mol.GetRingInfo().AtomRings()
    if len(cycles) == 0:
        maximum_ring_size = 0
    else:
        maximum_ring_size = max([len(ci) for ci in cycles])
    return maximum_ring_size



def passes_wehi_mcf(smi):
    mol =  Chem.MolFromSmiles(smi)
    h_mol = Chem.AddHs(mol)
    if any(h_mol.HasSubstructMatch(smarts) for smarts in _filters):
        return False 
    else: 
        return True
    
def pains_filt(mol):

    for k,v in dic.items():
        subs = Chem.MolFromSmarts( k )
        if subs != None:
            if mol.HasSubstructMatch( subs ):
                mol.SetProp(v,k)
    return [ prop for prop in mol.GetPropNames() ]



def substructure_violations(mol):
    """
    Check for substructure violates
    Return True: contains a substructure violation
    Return False: No substructure violation
    """
    violation = False

    forbidden_fragments = [
        "*1=**=*1",
        "*1*=*=*1",
        "*1~*=*1",
        "[F,Cl,Br]C=[O,S,N]",
        "[Br]-C-C=[O,S,N]",
        "[N,n,S,s,O,o]C[F,Cl,Br]",
        "[I]",
        "[S&X3]",
        "[S&X5]",
        "[S&X6]",
        "[B,N,n,O,S]~[F,Cl,Br,I]",
        "*=*=*=*",
        "*=[NH]",
        "[P,p]~[F,Cl,Br]",
        "SS", 
        "C#C", 
        "C=C=C", 
        "C=C=N", 
        "NNN", 
        "[*;R1]1~[*]~[*]~[*]1", 
        "OOO", 
        "[#8]1-[#6]2[#8][#6][#8][#6]12",                                                      # Epoxide group
        "N=C=O",                                                                              # Isocyanate group
        "C1CN1",                                                                              # Aziridine group 
        "[#6](=[#8])[F,Cl,Br,I]",                                                             # Acyl halides
        "[#6](=[#8])=[#6](-[#8])-[#6](=[#8])~[#8]",                                           # Quinone
        "N(-[#6])=[#7]-[#8]"                                                                  # Nitrosamine
    ]
    

    for ni in range(len(forbidden_fragments)):

        if mol.HasSubstructMatch(rdc.MolFromSmarts(forbidden_fragments[ni])) == True:
            # print('Violation frag: {} smi: {}'.format(forbidden_fragments[ni],  Chem.MolToSmiles(mol)) )
            violation = True
            break
        else:
            continue

    return violation


def filter_by_pattern(mol, pattern):
    """
    Check for presence of SMARTS pattern
    Return True: contains the pattern
    Return False: does not contain the pattern
    """
    violation = False

    if mol.HasSubstructMatch(rdc.MolFromSmarts(pattern)) == True:
        violation = True

    return violation

def filter_phosphorus(mol):
    """
    Check for presence of phopshorus fragment
    Return True: contains proper phosphorus
    Return False: contains improper phosphorus
    """
    violation = False

    if mol.HasSubstructMatch(rdc.MolFromSmarts("[P,p]")) == True:
        if mol.HasSubstructMatch(rdc.MolFromSmarts("*~[P,p](=O)~*")) == False:
            violation = True

    return violation
    
    
def apply_filters(smi):

    try: 
        if 'Si' in smi or 'Sn' in smi: # Atoms not appropriate for docking calculations 
            return False
        
        mol = smi2mol(smi)
        
        # Added after GDB-13 was filtered to get rid charged molecules
        if rdcmo.GetFormalCharge(mol) != 0:
            return False
        # Added after GDB-13 was filtered to get rid radicals
        elif rdcd.NumRadicalElectrons(mol) != 0:
            return False
        # Filter by bridgehead atoms
        elif rdcmd.CalcNumBridgeheadAtoms(mol) > 2:
            return False
        # Filter by ring size
        elif maximum_ring_size(mol) > 8:
            return False
        # Filter by proper phosphorus
        elif filter_phosphorus(mol):
            return False
        elif substructure_violations(mol):
            return False
        elif Descriptors.NumRotatableBonds(mol) > 10: 
            return False
        elif Descriptors.TPSA(mol) > 140: 
            return False 
        elif passes_wehi_mcf(smi) == False: 
            return False
        elif len(pains_filt(mol)) != 0: 
            return False
        else: 
            return True 
    except: 
        return False 
    
    


# Define a function to check if the SA score is < 4.5 and QED score is > 0.3

def is_good_molecule(smiles):
    # Convert the SMILES string to an RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False
    # Calculate the SA score
    sa_score = sascorer.calculateScore(mol)
    # print('SA score is: ', sa_score)
    # Calculate the QED score
    qed_score = QED.qed(mol)
    # Check if both conditions are met
    if sa_score < 4.5 and qed_score > 0.3:
        return True
    else:
        return False
    

def process_molecule(smi):
    if is_good_molecule(smi) and apply_filters(smi) and lipinski_filter(smi):
        return (smi, 'PASS')
    else:
        return (smi, 'FAIL')



if __name__ == '__main__':
    smi = 'C1=NC2=C(N=C1)C(=CC=N2)C1=CC=NC2=C1N=CN=C2'
    pass_filt = process_molecule(smi)
    print('Molecule: ', pass_filt)


