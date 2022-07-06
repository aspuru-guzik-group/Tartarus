from rdkit import Chem
from rdkit.Chem import AllChem

from rdkit import rdBase
rdBase.DisableLog('rdApp.error')

import numpy as np

import subprocess
import os

import time

from openbabel import pybel as pb
import pathlib as pl
cwd = pl.Path.cwd()
crest = pl.Path(os.environ["XTBHOME"] / pl.Path("bin") / pl.Path("crest"))
xtb = pl.Path(os.environ["XTBHOME"] / pl.Path("bin") / pl.Path("xtb"))

import subprocess as sp
import warnings
warnings.filterwarnings("ignore", message="Using default_file_mode other than 'r' is deprecated.")

from pyscf import gto, scf, dft, tddft
electronvolt = 27.211386245988

class molecule:
    def __init__(self, smiles, name, file=''):
        self.name = name
        self.smiles = smiles
        self.file = file
        self.generate_coordinates()
        return
    
    def generate_coordinates_rdkit(self):
        """
        convert smiles to xyz coordinates using rdkit
        Authors: Pascal Friederich, Robert Pollice
        """
        try:
            self.mol = Chem.MolFromSmiles(self.smiles)
        except:
            print("ERROR: could not convert %s to rdkit molecule."%(self.smiles))
            exit()
        try:
            self.mol = Chem.AddHs(self.mol)
        except:
            print("ERROR: could not add hydrogen to rdkit molecule of %s."%(self.smiles))
            exit()
        try:
            AllChem.EmbedMolecule(self.mol)
        except:
            print("ERROR: could not calculate 3D coordinates from rdkit molecule %s."%(self.smiles))
            exit()
        try:
            AllChem.MMFFOptimizeMolecule(self.mol)
        except:
            print("ERROR: could not optimize 3D coordinates for rdkit molecule %s."%(self.smiles))
            exit()
        try:
            Chem.rdMolTransforms.CanonicalizeMol(self.mol, normalizeCovar=True, ignoreHs=False)
        except:
            print("ERROR: could not canonicalize 3D coordinates for rdkit molecule %s."%(self.smiles))
            exit()
        try:
            block=Chem.MolToMolBlock(self.mol)
            blocklines=block.split("\n")
            self.xyz=[]
            self.atoms=[]
            for line in blocklines[4:]:
                if len(line.split())==4:
                    break
                self.atoms.append(line.split()[3])
                self.xyz.append([float(line.split()[0]),float(line.split()[1]),float(line.split()[2])])
            self.xyz=np.array(self.xyz)
            self.atom_count = len(self.atoms)
        except:
            print("ERROR: could not read xyz coordinates from rdkit molecule %s."%(self.smiles))
            exit()
        return
    
    def generate_coordinates_obabel(self):
        """
        convert smiles to xyz coordinates using openbabel
        Author: Robert Pollice
        """
        try:
            self.mol = pb.readstring('smi', self.smiles)
        except:
            print("ERROR: could not convert %s to openbabel molecule."%(self.smiles))
            exit()
        try:
            self.mol.make3D()
        except:
            print("ERROR: could not calculate 3D coordinates from openbabel molecule %s."%(self.smiles))
            exit()
        try:
            block = self.mol.write('xyz')
            blocklines=block.split("\n")
            self.xyz=[]
            self.atoms=[]
            for line in blocklines[2:-1]:
                self.atoms.append(line.split()[0])
                self.xyz.append([float(line.split()[1]),float(line.split()[2]),float(line.split()[3])])
            self.xyz=np.array(self.xyz)
            self.atom_count = len(self.atoms)
        except:
            print("ERROR: could not read xyz coordinates from openbabel molecule %s."%(self.smiles))
            exit()
        return
    
    def generate_coordinates_from_file(self):
        """
        read xyz coordinates from file
        """
        try:
            with open(pl.Path(cwd / self.file), 'r') as f:
                block = f.read()
                blocklines=block.split("\n")
                self.xyz=[]
                self.atoms=[]
                for line in blocklines[2:-1]:
                    self.atoms.append(line.split()[0])
                    self.xyz.append([float(line.split()[1]),float(line.split()[2]),float(line.split()[3])])
                self.xyz=np.array(self.xyz)
                self.atom_count = len(self.atoms)
        except:
            print("ERROR: could not generate xyz coordinates for molecule %s from file. Exiting."%(self.smiles))
            exit()
        return
    
    def generate_coordinates(self):
        # Read XYZ from file if specified
        if len(self.file) != 0:
            self.generate_coordinates_from_file()
        # Generate XYZ from SMILES if no file specified
        else:
            energies = []
            try:
                self.generate_coordinates_obabel()
                self.write_xyz()
                energies.append(self.xtb_energy(self.file))
                sp.run(['mv', self.file, 'obabel.xyz'], stdout=sp.DEVNULL, stderr=sp.DEVNULL)
            except:
                energies.append(10000.)
            try:
                self.generate_coordinates_rdkit()
                self.write_xyz()
                energies.append(self.xtb_energy(self.file))
                sp.run(['mv', self.file, 'rdkit.xyz'], stdout=sp.DEVNULL, stderr=sp.DEVNULL)
            except:
                energies.append(10000.)
            # both geometries failed
            if energies[0] == 10000. and energies[1] == 10000.:
                print("ERROR: could not generate xyz coordinates for molecule %s from file. Exiting."%(self.smiles))
                exit()
            # obabel coordinates are lower in energy
            elif np.argmin(energies) == 0:
                sp.run(['mv', 'obabel.xyz', self.file], stdout=sp.DEVNULL, stderr=sp.DEVNULL)
                sp.run(['rm', 'rdkit.xyz'], stdout=sp.DEVNULL, stderr=sp.DEVNULL)
            # rdkit coordinates are lower in energy
            elif np.argmin(energies) == 1:
                sp.run(['mv', 'rdkit.xyz', self.file], stdout=sp.DEVNULL, stderr=sp.DEVNULL)
                sp.run(['rm', 'obabel.xyz'], stdout=sp.DEVNULL, stderr=sp.DEVNULL)
            self.generate_coordinates_from_file()
            # clean up gfnff files
            sp.run(['rm', 'gfnff_charges', 'gfnff_topo'], stdout=sp.DEVNULL, stderr=sp.DEVNULL)
        return
    
    def xtb_energy(self, file):
        # Single point
        sp.run(['rm', 'gfnff_charges', 'gfnff_topo'], stdout=sp.DEVNULL, stderr=sp.DEVNULL)
        output = sp.run([xtb, file, '--gfnff'], stdout=sp.PIPE, stderr=sp.DEVNULL)
        # Read output
        line = [li for li in output.stdout.decode("utf-8").splitlines() if 'TOTAL ENERGY' in li]
        line = line[0].split()
        energy = float(line[3])
        return energy

    def write_xyz(self):
        text = str(self.atom_count)
        text += "\n"
        text += "Generated by tadf.py\n"
        for line in range(self.atom_count):
            text += self.atoms[line]
            text += format(format(self.xyz[line][0], '.9f'), '>14s')
            text += format(format(self.xyz[line][1], '.9f'), '>14s')
            text += format(format(self.xyz[line][2], '.9f'), '>14s')
            text += "\n"
        self.file = self.name + '.xyz'
        with open(pl.Path(cwd / self.file), 'w') as f:
            f.write(text)
        return


class computation:
    def __init__(self, file, name, atom_count):
        self.file = file
        self.name = name
        self.atom_count = atom_count
        # Clean topology files before running the workflow
        self.clean_files(thorough=False)
        self.pre_xtb()
        self.clean_files(thorough=False)
        self.crest()
        self.clean_files(thorough=False)
        self.xtb()
        self.clean_files(thorough=True)
        self.pyscf()
        os.system(f'rm -rf MRMSD {self.file} gfnff_adjacency')
        return
    
    def pre_xtb(self):
        # Initial geometry optimization needed for stable crest
        sp.run([xtb, self.file, '--gfn', '0', '--opt', 'normal', '--iterations', '4000'], stdout=sp.DEVNULL, stderr=sp.DEVNULL)
        return
        
    def crest(self):
        # Conformer search
        sp.run([crest, 'xtbopt.xyz', '-gff', '-mquick', '--noreftopo'], stdout=sp.DEVNULL, stderr=sp.DEVNULL)
        return
    
    def xtb(self):
        # Get best conformer from crest_best.xyz to avoid rare bug of crest where crest_best.xyz is nonsense
        with open('crest_best.xyz', 'w') as file:
            sp.run(['head', '-n', str(self.atom_count + 2), 'crest_conformers.xyz'], stdout=file, stderr=sp.DEVNULL)
        # Geometry optimization
        sp.run([xtb, 'crest_best.xyz', '--gfn', '0', '--opt', 'normal', '--iterations', '4000'], stdout=sp.DEVNULL, stderr=sp.DEVNULL)
        sp.run(['mv', 'xtbopt.xyz', self.file], stdout=sp.DEVNULL, stderr=sp.DEVNULL)
        return
    
    def clean_files(self, thorough=True):
        # Clean files
        if thorough == True:
            sp.run(['rm', 'bondlengths', 'charges', 'coord', 'coord.original', 'cregen_0.tmp', 'cregen_1.tmp', 'cre_members', 'crest_best.xyz', 'crest_conformers.xyz', 'crest.energies', 'crest_rotamers.xyz', 'gfnff_charges', 'gfnff_topo', '.history.xyz', 'struc.xyz', 'wbo', 'xtbopt.log', '.xtboptok', 'xtbrestart', 'xtbtopo.mol', 'xtblast.xyz'], stdout=sp.DEVNULL, stderr=sp.DEVNULL)
        elif thorough == False:
            sp.run(['rm', 'gfnff_charges', 'gfnff_topo'], stdout=sp.DEVNULL, stderr=sp.DEVNULL)
        return
    
    def pyscf(self):
        # TD-DFT single point
        self.mole = gto.Mole()
        self.mole.atom = self.file
        self.mole.basis = '6-31G*'
        self.mole.basis = {'H': '6-31G*', 'He': '6-31G*', 'Li': '6-31G*', 'Be': '6-31G*', 'B': '6-31G*', 'C': '6-31G*', 'N': '6-31G*', 'O': '6-31G*', 'F': '6-31G*', 'Ne': '6-31G*', 'Na': '6-31G*', 'Mg': '6-31G*', 'Al': '6-31G*', 'Si': '6-31G*', 'P': '6-31G*', 'S': '6-31G*', 'Cl': '6-31G*', 'Ar': '6-31G*', 'K': '6-31G*', 'Ca': '6-31G*', 'Sc': '6-31G*', 'Ti': '6-31G*', 'V': '6-31G*', 'Cr': '6-31G*', 'Mn': '6-31G*', 'Fe': '6-31G*', 'Co': '6-31G*', 'Ni': '6-31G*', 'Cu': '6-31G*', 'Zn': '6-31G*', 'Ga': '6-31G*', 'Ge': '6-31G*', 'As': '6-31G*', 'Se': '6-31G*', 'Br': '6-31G*', 'Kr': '6-31G*', 'I': 'lanl2dz'}
        self.mole.ecp = 'lanl2dz'
        self.mole.build()
        # Singlets
        self.method = dft.RKS(self.mole).density_fit(auxbasis='def2-universal-jkfit')
        self.method.xc = 'B3LYP'
        self.method.max_cycle = 512
        self.method.grids.level = 0
        self.method.conv_tol = 1E-7
        self.method.kernel()
        self.calculation = tddft.TDA(self.method)
        self.calculation.nstates = 2
        self.calculation.kernel()
        self.excitation_energies_singlets = self.calculation.e * electronvolt
        self.oscillator_strength_singlets = self.calculation.oscillator_strength(gauge='length')
        # Triplets
        self.calculation.singlet = False
        self.calculation.kernel()
        self.excitation_energies_triplets = self.calculation.e * electronvolt
        # Gather results
        self.fluorescence_energy = min(self.excitation_energies_singlets)
        self.singlet_triplet_gap = min(self.excitation_energies_singlets) - min(self.excitation_energies_triplets)
        self.oscillator_strength = self.oscillator_strength_singlets[0]
        self.results = np.array([self.fluorescence_energy, self.singlet_triplet_gap, self.oscillator_strength])
        return
    

def get_properties(smi): 
    
    start_time = time.time()
    try: 
        name_ = 'somet'
        mol = molecule(smi, 'random_name')
        mol.write_xyz()

        comp     = computation(mol.file, mol.name, mol.atom_count)
        results_ = comp.results # results_[0]=color; results_[1]=S-T; results_[2]=Oscillator Strength;        
        
        os.system('rm -rf core.*')
        run_time = time.time() - start_time
        
        if results_[1] <= 0.0 or results_[0]<=0.0 or results_[2]<=0.0: 
            return -10**4, -10**4, -10**4
        
        st, osc, combined = -results_[1], results_[2], results_[2]-results_[1]-np.abs(results_[0]-3.2)
        return st, osc, combined

    except: 
        return -10**4, -10**4, -10**4


if __name__ == '__main__': 
    smi = 'c1ccccc1'
    st, osc, combined= get_properties(smi)
    print(f'Singlet-triplet: {st}')
    print(f'Oscillator strength: {osc}')
    print(f'Combined obj: {combined}')


