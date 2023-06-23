# -*- coding: utf-8 -*-
"""
Created on Fri Jun  9 12:58:37 2023

@author: aksha
"""

import os 
import uuid
import time
import subprocess
import itertools
import argparse
import multiprocessing    

from filter_ import process_molecule 

def check_energy(lig_): 
    """
    Check the quality of a generated structure by computing its total energy using the Open Babel obenergy tool.
    Parameters:
        lig_ (str): the name of the ligand file in PDBQT format.
    Returns:
        total_energy (float): the computed total energy of the ligand in Kcal/mol.
    """
    # Check the quality of generated structure (some post-processing quality control):
    try: 
        ob_cmd = ['obenergy', lig_]
        command_obabel_check = subprocess.run(ob_cmd, capture_output=True)
        command_obabel_check = command_obabel_check.stdout.decode("utf-8").split('\n')[-2]
        total_energy         = float(command_obabel_check.split(' ')[-2])
    except: 
        total_energy = 10000 # Calculation has failed. 
        
    return total_energy


    

def run_docking_1syh(lig_location, out_location, method='qvina'): 
    """
    Perform molecular docking with a specific methods (QuickVina/Smina) on the 1SYH protein. 
    An exhaustiveness of 10 is used for the QuickVina calculations, while an 
    exhaustivesness of 100 is used for a smina calculation 

    Parameters
    ----------
    method : str, The calculation type to be run qvina/smina 

    Returns
    -------
    (float) Docking score.
    """
    if method == 'qvina': 
        command_run = subprocess.run(["./data/qvina", "--receptor", "./docking_structures/1syh/prot.pdbqt", "--ligand", lig_location, "--center_x", "21.492800140380858", "--center_y", "13.457733376820881", "--center_z", "23.175899950663247", "--size_x", "20", "--size_y", "20", "--size_z", "20", "--exhaustiveness", "10", "--out", out_location], capture_output=True)
    elif method == 'smina': 
        command_run = subprocess.run(["./data/smina", "--receptor", "./docking_structures/1syh/prot.pdbqt", "--ligand", lig_location, "--center_x", "21.492800140380858", "--center_y", "13.457733376820881", "--center_z", "23.175899950663247", "--size_x", "20", "--size_y", "20", "--size_z", "20", "--exhaustiveness", "100", "--out", out_location], capture_output=True)
    else: 
        raise Exception('Possible docking softwares: qvina/smina')

    # Ensure the pose of the output molecule is not broken: 
    pose_energy = check_energy(out_location)
    if pose_energy == 10000: # broken molecule generated (docking failed)
        return 10000
        
    # Obtain the docking score: 
    command_run = command_run.stdout.decode("utf-8").split('\n')

    docking_score = []
    for item in command_run: 
        line_split = item.split(' ')
        line_split = [x for x in line_split if x != '']
        if len(line_split) == 4: 
            try: 
                _ = float(line_split[0])
                vr_2 = float(line_split[1])
                _ = float(line_split[2])
                _ = float(line_split[3])
                docking_score.append(vr_2)
            except: 
                continue
    docking_score = min(docking_score)

    return docking_score


def run_docking_4lde(lig_location, out_location, method='qvina'): 
    """
    Perform molecular docking with a specific methods (QuickVina/Smina) on the 4LDE protein. 
    An exhaustiveness of 10 is used for the QuickVina calculations, while an 
    exhaustivesness of 100 is used for a smina calculation 

    Parameters
    ----------
    method : str, The calculation type to be run qvina/smina 

    Returns
    -------
    (float) Docking score.
    """
    if method == 'qvina': 
        command_run = subprocess.run(["./data/qvina", "--receptor", "./docking_structures/4lde/prot.pdbqt", "--ligand", lig_location, "--center_x", "-2.942962976793448", "--center_y", "-12.915592617458767", "--center_z", "-50.99233344749168", "--size_x", "20", "--size_y", "20", "--size_z", "20", "--exhaustiveness", "10", "--out", out_location], capture_output=True)
    elif method == 'smina': 
        command_run = subprocess.run(["./data/smina", "--receptor", "./docking_structures/4lde/prot.pdbqt", "--ligand", lig_location, "--center_x", "-2.942962976793448", "--center_y", "-12.915592617458767", "--center_z", "-50.99233344749168", "--size_x", "20", "--size_y", "20", "--size_z", "20", "--exhaustiveness", "100", "--out", out_location], capture_output=True)
    else: 
        raise Exception('Possible docking softwares: qvina/smina')

    # Ensure the pose of the output molecule is not broken: 
    pose_energy = check_energy(out_location)
    if pose_energy == 10000: # broken molecule generated (docking failed)
        return 10000
        
    # Obtain the docking score: 
    command_run = command_run.stdout.decode("utf-8").split('\n')

    docking_score = []
    for item in command_run: 
        line_split = item.split(' ')
        line_split = [x for x in line_split if x != '']
        if len(line_split) == 4: 
            try: 
                _ = float(line_split[0])
                vr_2 = float(line_split[1])
                _ = float(line_split[2])
                _ = float(line_split[3])
                docking_score.append(vr_2)
            except: 
                continue
    docking_score = min(docking_score)

    return docking_score


def run_docking_6y2f(lig_location, out_location, method='qvina'): 
    """
    Perform molecular docking with a specific methods (QuickVina/Smina) on the 6Y2F protein. 
    An exhaustiveness of 10 is used for the QuickVina calculations, while an 
    exhaustivesness of 100 is used for a smina calculation 

    Parameters
    ----------
    method : str, The calculation type to be run qvina/smina 

    Returns
    -------
    (float) Docking score.
    """
    if method == 'qvina': 
        command_run = subprocess.run(["./datasets/qvina", "--receptor", "./docking_structures/6y2f/prot.pdbqt", "--ligand", lig_location, "--center_x", "11.026168696851615", "--center_y", "-0.6082891440804464", "--center_z", "20.840999947973046", "--size_x", "10", "--size_y", "10", "--size_z", "10", "--exhaustiveness", "10", "--out", out_location], capture_output=True)
    elif method == 'smina': 
        command_run = subprocess.run(["./datasets/smina", "--receptor", "./docking_structures/6y2f/prot.pdbqt", "--ligand", lig_location, "--center_x", "11.026168696851615", "--center_y", "-0.6082891440804464", "--center_z", "20.840999947973046", "--size_x", "10", "--size_y", "10", "--size_z", "10", "--exhaustiveness", "100", "--out", out_location], capture_output=True)
    else: 
        raise Exception('Possible docking softwares: qvina/smina')

    # Ensure the pose of the output molecule is not broken: 
    pose_energy = check_energy(out_location)
    if pose_energy == 10000: # broken molecule generated (docking failed)
        return 10000
        
    # Obtain the docking score: 
    command_run = command_run.stdout.decode("utf-8").split('\n')

    docking_score = []
    for item in command_run: 
        line_split = item.split(' ')
        line_split = [x for x in line_split if x != '']
        if len(line_split) == 4: 
            try: 
                _ = float(line_split[0])
                vr_2 = float(line_split[1])
                _ = float(line_split[2])
                _ = float(line_split[3])
                docking_score.append(vr_2)
            except: 
                continue
    docking_score = min(docking_score)

    return docking_score


def generate_unique_file_name(base_name, extension):
    timestamp = int(time.time() * 1000)
    unique_id = uuid.uuid4().hex
    file_name = f"{base_name}_{timestamp}_{unique_id}.{extension}"
    return file_name


def perform_calc_single(smi, receptor_type, docking_program='qvina'): 
    """
    Performs docking calculations on a single molecule-receptor pair.
    
    This function takes a molecule (in SMILES format), a receptor type, and an optional docking program 
    (default: 'qvina' / 'smina'). It generates a 3D structure for the molecule, performs energy check, and 
    runs a docking simulation with the receptor specified. It returns the docking score.
    
    Parameters
    ----------
    smi : str
        The SMILES string of the molecule.
    
    receptor_type : str
        The type of the receptor. This function currently supports '1syh', '4lde', and '6y2f'.
    
    docking_program : str, optional
        The docking software to be used for the docking calculation. Default is 'qvina'.
    
    Returns
    -------
    float
        The docking score of the molecule with the receptor. If the process fails or the 
        molecule's energy is too high, the function returns a value of 10^4.
    
    Raises
    ------
    It does not raise any exceptions, but rather catches them and returns a docking score of 10^4.
    
    Notes
    -----
    The function employs Open Babel ('obabel') for the 3D structure generation of the molecule 
    and the specified docking program for the docking process. Intermediate files (pdbqt format) 
    created during the execution are deleted before the function returns.
    """
    try: 
    
        pass_filt =  process_molecule(smi)
        if pass_filt[1] == 'Fail': 
            return 10**4
    
        output_filename = generate_unique_file_name('lig', 'pdbqt')
        # print('smi: {} fname: {}'.format(smi, output_filename))
        cmd = ["obabel", "-ismi","-:" + smi,"-O", output_filename, "--gen3d"]
        subprocess.run(cmd, timeout=20)
        
        # Ensure a stable molecule: 
        lig_energy = check_energy(output_filename)
    
        # Specifying docking input file & output file: 
        lig_location = output_filename
        out_location = generate_unique_file_name('pose', 'pdbqt')
        
        # Perform docking: 
        if lig_energy < 10000: 
            if receptor_type == '1syh': 
                score_ = run_docking_1syh(lig_location, out_location, method=docking_program)
            if receptor_type == '4lde': 
                score_ = run_docking_4lde(lig_location, out_location, method=docking_program)
            if receptor_type == '6y2f': 
                score_ = run_docking_6y2f(lig_location, out_location, method=docking_program)
        else: 
            return 10**4
    
                
        os.system('rm {} {}'.format(output_filename, out_location))
        
    except: 
        os.system('rm {} {}'.format(output_filename, out_location))
        return 10**4

    return score_


if __name__ == '__main__': 
    
    smi = 'C1=NC2=C(N=C1)C(=CC=N2)C1=CC=NC2=C1N=CN=C2'
    
    # Run calculation for the 1syh recepotor: 
    score = perform_calc_single(smi, '1syh', docking_program='qvina')

    # Run calculation for the 4lde recepotor: 
    score = perform_calc_single(smi, '4lde', docking_program='qvina')

    # Run calculation for the 6y2f recepotor: 
    score = perform_calc_single(smi, '6y2f', docking_program='qvina')




