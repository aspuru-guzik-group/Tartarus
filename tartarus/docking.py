import os 
import time 
import multiprocessing
import inspect
import tempfile
from pathlib import Path
from .utils import run_command

from rdkit import Chem 
from rdkit.Chem import MolFromSmiles
from rdkit.Chem.Crippen import MolLogP
from rdkit.Chem.Descriptors import ExactMolWt
from rdkit.Chem.Lipinski import NumHAcceptors, NumHDonors
# from filter_ import maximum_ring_size, filter_phosphorus, substructure_violations

from rdkit.Chem import MolFromSmiles as smi2mol
import rdkit.Chem.rdmolops as rdcmo
import rdkit.Chem.Descriptors as rdcd
import rdkit.Chem.rdMolDescriptors as rdcmd

import rdkit.Chem as rdc


def lipinski_filter(smiles):
    mol = MolFromSmiles(smiles)
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


def substructure_violations(mol):
    """
    Check for substructure violates
    Return True: contains a substructure violation
    Return False: No substructure violation
    """
    violation = False
    # forbidden_fragments = ['*1=**=*1', '*1*=*=*1', '*1~*=*1', '[F,Cl,Br]C=[O,S,N]','[Cl,Br]-C-C=[O,S,N]','[N,n,S,s,O,o]C[F,Cl,Br]','[I]','*=[S,s;!R]', '[S&X3]', '[S&X4]', '[S&X5]', '[S&X6]', '[N+]', '[B,N,n,O,S]~[F,Cl,Br,I]', '*=*=*', '*#*', '[O,o,S,s]~[O,o,S,s]', '[N,n,O,o,S,s]~[N,n,O,o,S,s]~[N,n,O,o,S,s]', '[C,c]~N=,:[O,o,S,s;!R]', '[N,n,O,o,S,s]~[N,n,O,o,S,s]~[C,c]=,:[O,o,S,s,N,n;!R]', '*=[NH]', '*=N-[*;!R]', '*~[N,n,O,o,S,s]-[N,n,O,o,S,s;!R]']
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
        "C=C=C"
    ]

    for ni in range(len(forbidden_fragments)):

        if mol.HasSubstructMatch(rdc.MolFromSmarts(forbidden_fragments[ni])) == True:
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
        else: 
            return True 
    except: 
        return False 


def get_score(smi, docking_target='1syh', verbose=False, scratch='/tmp'): 
    # Create and switch to temporary directory
    owd = Path.cwd()
    scratch_path = Path(scratch)
    tmp_dir = tempfile.TemporaryDirectory(dir=scratch_path)
    os.chdir(tmp_dir.name)

    system = lambda x: run_command(x, verbose)
    try: 
        with open('test.smi', 'w') as f: 
            f.writelines(smi)
        system('obabel test.smi --gen3D -O hit1.pdb')
        
        cpu_count = multiprocessing.cpu_count()

        path = os.path.join(
            os.path.dirname(inspect.getfile(apply_filters)), 
            'docking_structures'
        )

        command = (
            f'{path}/smina.static -r {path}/{docking_target}/prot.pdb '
            f'-l hit1.pdb --autobox_ligand {path}/{docking_target}/lig_.pdb '
            '--autobox_add 2 --scoring vina -o dock_1.pdb --log dock_1.log '
            f'--cpu {cpu_count} --num_modes 1 --exhaustiveness 10'
        )
        system(command)
        
        # Read in the results:  
        with open('./dock_1.log', 'r') as f: 
            lines = f.readlines()
        lines = lines[25: ]
        
        scores = []
        for line in lines: 
            A = line.split(' ')
            A = [x for x in A if x != '']
            scores.append(float(A[1]))
    except: 
        scores = 10**4

    # Remove temporary directory
    os.chdir(owd)
    tmp_dir.cleanup()
    
        
    # Delete log files: 
    system('rm dock_1.log dock_1.pdb test.smi hit1.pdb')
    
    return min(scores)


def get_1syh_score(smi: str, verbose: bool = False): 
    if apply_filters(smi) == True and lipinski_filter( smi ) == True: 
        score = get_score(smi=smi, docking_target='1syh', verbose=verbose) # '4lde', '6y2f', '1syh'
        return score # MINUS SIGN FOR OPTIMIZATION! 
    else: 
        return 10**4

def get_6y2f_score(smi: str, verbose: bool = False): 
    if apply_filters(smi) == True and lipinski_filter( smi ) == True: 
        score = get_score(smi=smi, docking_target='6y2f', verbose=verbose) # '4lde', '6y2f', '1syh'
        return score # MINUS SIGN FOR OPTIMIZATION! 
    else: 
        return 10**4

def get_4lde_score(smi: str, verbose: bool = False): 
    if apply_filters(smi) == True and lipinski_filter( smi ) == True: 
        score = get_score(smi=smi, docking_target='4lde', verbose=verbose) # '4lde', '6y2f', '1syh'
        return score # MINUS SIGN FOR OPTIMIZATION! 
    else: 
        return 10**4

if __name__ == '__main__':        
    smi = 'OC1=CC2=C(C=CC=C2)C(=C1)C1=CC(O)=CC2=C1C=CN=C2'
    score = get_4lde_score(smi)
    print(f'4lde score: {score}')

    score = get_6y2f_score(smi)
    print(f'6y2f score: {score}')

    score = get_1syh_score(smi)
    print(f'1syh score: {score}')

