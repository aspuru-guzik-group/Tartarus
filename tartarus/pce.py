import os, sys
import numpy as np
import inspect

import rdkit
from rdkit.Chem import AllChem, DataStructs
from rdkit.Chem import AllChem as Chem

from rdkit.Chem import RDConfig
sys.path.append(os.path.join(RDConfig.RDContribDir, 'SA_Score'))
import sascorer

import torch
import torch.nn as nn


# Surrogate model for Jsc
def gaussian(x, A, B):
    return A * np.exp(-x** 2 / B)

def get_properties(smile): 
    # Create mol object
    mol = Chem.MolFromSmiles(smile)
    mol = Chem.AddHs(mol)
    if mol == None: 
        return "INVALID"
    charge = Chem.rdmolops.GetFormalCharge(mol)
    atom_number = mol.GetNumAtoms()

    sas = sascorer.calculateScore(mol)
    
    with open('test.smi', 'w') as f: 
        f.writelines([smile])

    # Prepare the input file: 
    os.system('obabel test.smi --gen3D -O test.xyz')

    # Run the preliminary xtb: 
    command_pre = 'CHARGE={};xtb {} --gfn 0 --opt normal -c $CHARGE --iterations 4000'.format(charge, 'test.xyz')
    os.system(command_pre)
    os.system("rm ./gfnff_charges ./gfnff_topo")

    # Run crest conformer ensemble
    command_crest = 'CHARGE={};crest {} -gff -mquick -chrg $CHARGE --noreftopo'.format(charge, 'xtbopt.xyz')
    os.system(command_crest)
    os.system('rm ./gfnff_charges ./gfnff_topo')
    os.system('head -n {} crest_conformers.xyz > crest_best.xyz'.format(atom_number+2))

    # Run the calculation: 
    command = 'CHARGE={};xtb {} --opt normal -c $CHARGE --iterations 4000 > out_dump'.format(charge, 'crest_best.xyz')
    os.system(command)

    # Read the output: 
    with open('./out_dump', 'r') as f: 
        text_content = f.readlines()

    output_index = [i for i in range(len(text_content)) if 'Property Printout' in text_content[i]]
    text_content = text_content[output_index[0]: ]
    homo_data = [x for x in text_content if '(HOMO)' in x]
    lumo_data = [x for x in text_content if '(LUMO)' in x]
    homo_lumo_gap = [x for x in text_content if 'HOMO-LUMO GAP' in x]
    mol_dipole    = [text_content[i:i+4] for i,x in enumerate(text_content) if 'molecular dipole:' in x]
    lumo_val      = float(lumo_data[0].split(' ')[-2])
    homo_val = float(homo_data[0].split(' ')[-2])
    homo_lumo_val  = float(homo_lumo_gap[0].split(' ')[-5])
    mol_dipole_val = float(mol_dipole[0][-1].split(' ')[-1])

    # Determine value of custom function for optimization
    HL_range_rest = homo_lumo_val # Good range for the HL gap: 0.8856-3.2627 
    if 0.8856 <= HL_range_rest <= 3.2627: 
        HL_range_rest = 1.0
    elif HL_range_rest < 0.8856: 
        HL_range_rest = 0.1144 + homo_lumo_val
    else: 
        HL_range_rest = 4.2627 - HL_range_rest
    function_ = mol_dipole_val + HL_range_rest - lumo_val # Maximize this function
    
    # Compute calibrated homo and lumo levels
    homo_cal = homo_val * 0.8051030400316004 + 2.5376777453204133
    lumo_cal = lumo_val * 0.8787863933542347 + 3.7912767464357200

    # Define parameters for Scharber model
    A = 433.11633173034136
    B = 2.3353220382662894
    Pin = 900.1393292842149

    # Scharber model objective 1: Optimization of donor for phenyl-C61-butyric acid methyl ester (PCBM) acceptor
    voc_1 = (abs(homo_cal) - abs(-4.3)) - 0.3
    if voc_1 < 0.0:
        voc_1 = 0.0
    lumo_offset_1 = lumo_cal + 4.3
    if lumo_offset_1 < 0.3:
        pce_1 = 0.0
    else:
        jsc_1 = gaussian(lumo_cal - homo_cal, A, B)
        if jsc_1 > 415.22529811760637:
            jsc_1 = 415.22529811760637
        pce_1 = 100 * voc_1 * 0.65 * jsc_1 / Pin

    # Scharber model objective 2: Optimization of acceptor for poly[N-90-heptadecanyl-2,7-carbazole-alt-5,5-(40,70-di-2-thienyl-20,10,30-benzothiadiazole)] (PCDTBT) donor
    voc_2 = (abs(-5.5) - abs(lumo_cal)) - 0.3
    if voc_2 < 0.0:
        voc_2 = 0.0
    lumo_offset_2 = -3.6 - lumo_cal
    if lumo_offset_2 < 0.3:
        pce_2 = 0.0
    else:
        jsc_2 = gaussian(lumo_cal - homo_cal, A, B)
        if jsc_2 > 415.22529811760637:
            jsc_2 = 415.22529811760637
        pce_2 = 100 * voc_2 * 0.65 * jsc_2 / Pin

    # Delete all the output files: 
    os.system('rm xtbopt.log xtbopt.xyz xtbrestart xtbtopo.mol charges out_dump test.smi test.xyz wbo .xtboptok')
    os.system('rm bondlengths coord coord.original cregen_0.tmp  cre_members crest_conformers.xyz crest.energies crest_rotamers.xyz struc.xyz .CHRG .history.xyz crest_best.xyz')
    os.system('rm -rf MRMSD gfnff_adjacency')

    return mol_dipole_val, homo_lumo_val, lumo_val, function_, pce_1, pce_2, sas




def get_fingerprint(smile, nBits, ecfp_degree=2):
    m1 = Chem.MolFromSmiles(smile)
    fp = AllChem.GetMorganFingerprintAsBitVect(m1,ecfp_degree, nBits=nBits)
    x = np.zeros((0, ), dtype=np.int8)
    DataStructs.ConvertToNumpyArray(fp, x)

    return x


class SurrogateModel(object):
    def __init__(self, model_list_dir, use_ensemble=False):
        super(SurrogateModel, self).__init__()
        self.use_ensemble = use_ensemble
        
        model_state_dicts = os.listdir(model_list_dir)
        self.model_list = []
        for model_state_dict in model_state_dicts:
            self.model_list.append(torch.load(os.path.join(model_list_dir, model_state_dict)))

            # print(torch.load(os.path.join(model_list_dir, model_state_dict)))
            torch.load(os.path.join(model_list_dir, model_state_dict))

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

        
def get_surrogate_properties(smiles):
    ''' Load surrogate models trained on dataset. Make predictions of material
    properties. 
    '''
    path = os.path.join(os.path.dirname(inspect.getfile(get_fingerprint)), 'pce_pretrained_models')
    try: 
        surrogate_model_homo_lumo = SurrogateModel(f'{path}/homo_lumo')
        surrogate_model_lumo_val = SurrogateModel(f'{path}/lumo_val')
        surrogate_model_dipole = SurrogateModel(f'{path}/dipole')
        surrogate_model_function = SurrogateModel(f'{path}/function')
        
        homo_lumo = surrogate_model_homo_lumo.forward(smiles) 
        lumo_val = surrogate_model_lumo_val.forward(smiles)
        dipole = surrogate_model_dipole.forward(smiles)
        function = surrogate_model_function.forward(smiles)
    
        return float(dipole), float(homo_lumo), float(-lumo_val), float(function)
    
    except: 
        return float(-10**4), float(-10**4), float(-10**4), float(-10**4)

if __name__ == '__main__':
    smi = 'c1ccccc1'
    dipole, hl_gap, lumo, obj, pce_1, pce_2, sas = get_properties(smi)
    print(f'Dipole: {dipole}')
    print(f'HOMO-LUMO Gap: {hl_gap}')
    print(f'LUMO: {lumo}')
    print(f'Combined obj: {obj}')
    print(f'PCE1: {pce_1}')
    print(f'PCE2: {pce_2}')
    print(f'SAS: {sas}')

    smi = 'c1ccccc1'
    dipole, hl_gap, lumo, obj = get_surrogate_properties(smi)
    print(f'Dipole: {dipole}')
    print(f'HOMO-LUMO Gap: {hl_gap}')
    print(f'LUMO: {lumo}')
    print(f'Combined obj: {obj}')


