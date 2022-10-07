"""Tools for reactivity benchmark."""

from __future__ import annotations

import os
from pathlib import Path
import subprocess
from subprocess import TimeoutExpired
import sys
import tempfile
import time

from morfeus.conformer import _add_conformers_to_mol
import numpy as np
from numpy.typing import NDArray
import polanyi.config
from polanyi.data import HARTREE_TO_KCAL
from polanyi.interpolation import interpolate_geodesic
from polanyi.io import read_xyz
from polanyi.workflow import crest_constrained, opt_ts_python, opt_xtb
from polanyi.xtb import opt_crest

import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import RDConfig
from scipy.spatial import distance_matrix
from sklearn.covariance import EllipticEnvelope

# Import SAScore
sys.path.append(os.path.join(RDConfig.RDContribDir, 'SA_Score'))
import sascorer

KEYWORDS_CREST = ["-mquick", "-gfn2//gfnff", "-noreftopo"]

def geometry_from_openbabel(mol: Chem.Mol) -> tuple[list, NDArray[float]]:
    """Get geometry from openbabel."""
    temp_file_mol = tempfile.NamedTemporaryFile()
    temp_file_xyz = tempfile.NamedTemporaryFile()
    Chem.MolToMolFile(mol, temp_file_mol.name)
    cmd = f"obabel -imol --gen3d -oxyz {temp_file_mol.name} -O {temp_file_xyz.name}.xyz"
    subprocess.run(cmd.split(), timeout=20, capture_output=True)
    elements, coordinates = read_xyz(f"{temp_file_xyz.name}.xyz")
    temp_file_xyz.close()
    temp_file_mol.close()
    return elements, coordinates


def apply_chiral_template(mol: Chem.Mol, pattern: Chem.Mol) -> Chem.Mol:
    """Apply chiral template to Mol object."""
    # Apply SMARTS pattern to Mol
    match = mol.GetSubstructMatch(pattern)

    # Find indices of matched bonds in Mol
    indices_matched = []
    for bond in pattern.GetBonds():
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        bond = mol.GetBondBetweenAtoms(match[i], match[j])
        indices_matched.append(bond.GetIdx())

    # Reorder bond indices to match SMARTS pattern
    indices_all = list(range(mol.GetNumBonds()))
    indices_rest = [i for i in indices_all if i not in indices_matched]
    indices_new = indices_matched + indices_rest

    # Create new Mol. Delete and re-add bonds in new order
    rw_mol = Chem.RWMol(mol)

    bond_info = []
    for bond in list(rw_mol.GetBonds()):
        bond_info.append(
            (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), bond.GetBondType())
        )
        rw_mol.RemoveBond(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())

    for i in indices_new:
        rw_mol.AddBond(*bond_info[i])

    # Set chiral tags from template to new Mol
    for i, atom in enumerate(pattern.GetAtoms()):
        chiral_tag = atom.GetChiralTag()
        rw_mol.GetAtomWithIdx(match[i]).SetChiralTag(chiral_tag)

    # Recover new Mol
    new_mol = rw_mol.GetMol()
    Chem.SanitizeMol(new_mol)

    return new_mol


def process_smiles(
    smiles: str,
) -> tuple[
    tuple[Chem.Mol, Chem.Mol, Chem.Mol], tuple[int, int, int], tuple[int, int, int]
]:
    """Process smiles to create reactant, product and TS with indices."""
    # Create mol
    mol = Chem.MolFromSmiles(smiles)

    # Create 3D strcture
    mol = Chem.AddHs(mol)
    try:
        _, coordinates = geometry_from_openbabel(mol)
        _add_conformers_to_mol(mol, [coordinates])
    except TimeoutExpired:
        return_code = AllChem.EmbedMolecule(
            mol,
            maxAttempts=mol.GetNumAtoms() * 1000,
            useSmallRingTorsions=True,
            useMacrocycleTorsions=True,
            useRandomCoords=True,
            forceTol=0.01,
            ETversion=2,
        )
        if return_code != 0:
            raise RuntimeError(f"Geometry embedding failed with return code: {return_code}")
    
    AllChem.MMFFOptimizeMolecule(mol)

    # Find matches for reactive core
    smarts = "C12CCC([*]2)C34C([*]5)C=CC5C13[*]4"
    pattern = Chem.MolFromSmarts(smarts)
    matches = mol.GetSubstructMatches(pattern)

    if len(matches) == 0:
        raise ValueError("No matches for reactive core.")
    elif len(matches) > 1:
        raise ValueError("More than one matches for reactive core.")
    match = matches[0]

    # Do conformer search
    elements = [atom.GetSymbol() for atom in mol.GetAtoms()]
    coordinates = mol.GetConformer().GetPositions()
    conformer_ensemble = opt_crest(elements, coordinates, keywords=KEYWORDS_CREST)
    mol.RemoveAllConformers()
    _add_conformers_to_mol(mol, [conformer_ensemble[0].coordinates])

    # Take out reactive carbon atoms
    index_c_u_l = match[9]
    index_c_u_r = match[8]
    index_c_s_l = match[1]
    index_c_s_r = match[2]
    atom_c_s_l = mol.GetAtomWithIdx(index_c_s_l)
    atom_c_s_r = mol.GetAtomWithIdx(index_c_s_r)

    # Find hydrogen atoms on unsaturade bridge based on distance matrix
    indices_h_l = [
        neighbor.GetIdx()
        for neighbor in atom_c_s_l.GetNeighbors()
        if neighbor.GetSymbol() == "H"
    ]
    indices_h_r = [
        neighbor.GetIdx()
        for neighbor in atom_c_s_r.GetNeighbors()
        if neighbor.GetSymbol() == "H"
    ]
    distance_matrix = Chem.Get3DDistanceMatrix(mol)
    index_h_l = indices_h_l[np.argmin(distance_matrix[indices_h_l, index_c_u_l])]
    index_h_r = indices_h_r[np.argmin(distance_matrix[indices_h_r, index_c_u_r])]

    # Create editable mol
    rw_mol = Chem.RWMol(mol)

    # Set bond orders
    bond = rw_mol.GetBondBetweenAtoms(index_c_u_l, index_c_u_r)
    bond.SetBondType(Chem.BondType.SINGLE)
    bond = rw_mol.GetBondBetweenAtoms(index_c_s_l, index_c_s_r)
    bond.SetBondType(Chem.BondType.DOUBLE)

    # Remove bonds between hydrogen atoms on one side and add them on the other side
    rw_mol.RemoveBond(index_c_s_l, index_h_l)
    rw_mol.RemoveBond(index_c_s_r, index_h_r)
    rw_mol.AddBond(index_c_u_l, index_h_l, Chem.BondType.SINGLE)
    rw_mol.AddBond(index_c_u_r, index_h_r, Chem.BondType.SINGLE)

    # Generate product mol in 3D
    mol_prod = rw_mol.GetMol()
    Chem.SanitizeMol(mol_prod)

    n_steps = 10
    h_dist = 1.1
    for dist_1, dist_2 in zip(
        np.linspace(distance_matrix[index_c_u_l, index_h_l] - h_dist, h_dist, n_steps),
        np.linspace(distance_matrix[index_c_u_r, index_h_r] - h_dist, h_dist, n_steps),
    ):
        ff_properties = AllChem.MMFFGetMoleculeProperties(mol_prod)
        ff = AllChem.MMFFGetMoleculeForceField(mol_prod, ff_properties)
        ff.AddDistanceConstraint(index_c_u_l, index_h_l, dist_1, dist_1, 10)
        ff.AddDistanceConstraint(index_c_u_r, index_h_r, dist_2, dist_2, 10)
        ff.Initialize()
        ff.Minimize()
    AllChem.MMFFOptimizeMolecule(mol_prod)

    mol_prod_ts = Chem.Mol(mol_prod)

    # Do conformer search
    elements = [atom.GetSymbol() for atom in mol_prod.GetAtoms()]
    coordinates = mol_prod.GetConformer().GetPositions()
    conformer_ensemble = opt_crest(elements, coordinates, keywords=KEYWORDS_CREST)
    mol_prod.RemoveAllConformers()
    _add_conformers_to_mol(mol_prod, [conformer_ensemble[0].coordinates])

    return (
        (mol, mol_prod_ts, mol_prod),
        (index_c_s_l, index_c_u_l, index_h_l),
        (index_c_s_r, index_c_u_r, index_h_r),
    )


def barrier_from_smiles(
    smiles: str, n_images: int = 5, i_image: int = None
) -> tuple[float, float]:
    """Calculate activation and reaction energy from SMILES."""
    # Generate geomertries of reactant and products
    (
        (mol_r, mol_p_ts, mol_p),
        (index_c_s_l, index_c_u_l, index_h_l),
        (index_c_s_r, index_c_u_r, index_h_r),
    ) = process_smiles(smiles)
    elements = [atom.GetSymbol() for atom in mol_r.GetAtoms()]
    coordinates_r = mol_r.GetConformer().GetPositions()
    coordinates_p = mol_p.GetConformer().GetPositions()
    coordinates_p_ts = mol_p_ts.GetConformer().GetPositions()

    # Optimize with XTB
    keywords = ["--opt", "--gfnff"]
    coordinates_r = opt_xtb(elements, coordinates_r, keywords=keywords)
    coordinates_p = opt_xtb(elements, coordinates_p, keywords=keywords)
    coordinates_p_ts = opt_xtb(elements, coordinates_p_ts, keywords=keywords)

    # Generate guess for transition state
    path = interpolate_geodesic(
        elements, [coordinates_r, coordinates_p_ts], n_images=n_images
    )
    if i_image is None:
        i_image = n_images // 2
    coordinates_guess = path[i_image]

    # Optimize TS
    maxsteps = 200
    conv_params = {  # These are the default settings
        "convergence_energy": 1e-3,  # Eh
        "convergence_grms": 3e-3,  # Eh/Bohr
        "convergence_gmax": 4.5e-3,  # Eh/Bohr
        "convergence_drms": 1.2e-3,  # Angstrom
        "convergence_dmax": 1.8e-3,  # Angstrom
        "coordsys": "cart",
    }

    kw_opt = {"conv_params": conv_params, "maxsteps": maxsteps}
    results = opt_ts_python(
        elements, [coordinates_r, coordinates_p], coordinates_guess, kw_opt=kw_opt
    )

    # Check convergence
    n_steps = len(results.coordinates_opt)
    if n_steps == maxsteps:
        raise RuntimeError(f"Optimization did not converge in {maxsteps} steps.")

    # Do conformer search and optimize TS again
    dm = distance_matrix(results.coordinates_opt, results.coordinates_opt)
    d_1 = dm[index_c_s_l, index_h_l]
    d_2 = dm[index_c_u_l, index_h_l]
    d_3 = dm[index_c_s_r, index_h_r]
    d_4 = dm[index_c_u_r, index_h_r]

    distance_constraints = {
        (index_c_s_l + 1, index_h_l + 1): d_1,
        (index_c_u_l + 1, index_h_l + 1): d_2,
        (index_c_s_r + 1, index_h_r + 1): d_3,
        (index_c_u_r + 1, index_h_r + 1): d_4,
    }

    conformer_ensemble_ts = crest_constrained(
        elements,
        results.coordinates_opt,
        keywords=KEYWORDS_CREST,
        distance_constraints=distance_constraints,
        fc=1.0,
    )
    coordinates_ts_crest = conformer_ensemble_ts[0].coordinates

    results = opt_ts_python(
        elements, [coordinates_r, coordinates_p], coordinates_ts_crest, kw_opt=kw_opt
    )

    # Check convergence
    n_steps = len(results.coordinates_opt)
    if n_steps == maxsteps:
        raise RuntimeError(f"Optimization did not converge in {maxsteps} steps.")

    # Calculate reaction and activation energy
    energy_reaction = results.shift_results.energy_diff_gfn * HARTREE_TO_KCAL

    energy_reactant = results.shift_results.energies_ff[0]
    energy_ts = min(results.opt_results.energies_diabatic[-1])
    energy_activation = (energy_ts - energy_reactant) * HARTREE_TO_KCAL

    # Generate SMILES with stereochemistry
    mol_stereo = Chem.Mol(mol_r)
    Chem.AssignStereochemistryFrom3D(mol_stereo)
    mol_stereo = Chem.RemoveHs(mol_stereo)
    smiles_stereo = Chem.MolToSmiles(mol_stereo)

    return energy_activation, energy_reaction, smiles_stereo


def run_reaction(
    smiles: str, n_procs: int = 1, scratch: str = "/tmp"
) -> tuple[str, float, float, float, float]:
    """Run one reaction."""
    # Create and switch to temporary directory
    owd = Path.cwd()
    scratch_path = Path(scratch)
    tmp_dir = tempfile.TemporaryDirectory(dir=scratch_path)
    os.chdir(tmp_dir.name)

    # Start timer
    start_time = time.time()

    # Set up polanyi
    polanyi.config.OMP_NUM_THREADS = str(n_procs)
    polanyi.config.TMP_DIR = scratch

    # Applying the right stereochemistry to the core
    smarts = "[H][C@@]1(*)[C@;R2](*)2[C@@]34[C@@;R2]5(*)[C;R1](*)=[C;R1](*)[C@;R2](*)([*;R2]5)[C@@]3([*;R1]4)[C@](*)([*;R2]2)[C@@;R1]1([*])[H]"
    pattern = Chem.MolFromSmarts(smarts)

    mol = Chem.MolFromSmiles(smiles)
    # Calculate SAScore
    sa_score = sascorer.calculateScore(mol)
    mol = Chem.AddHs(mol)
    mol = apply_chiral_template(mol, pattern)
    mol = Chem.RemoveHs(mol)
    smiles = Chem.MolToSmiles(mol)

    # Run the calculation
    failed = False
    try: 
        activation_energy, reaction_energy, smiles_stereo = barrier_from_smiles(smiles)
    except Exception:
        print('Reactivity simulation has failed for this molecule.')
        failed = True
        smiles_stereo = None
        activation_energy = 10**4
        reaction_energy = 10**4
        sa_score = 10**4

    # Calculate run time
    end_time = time.time()
    run_time = end_time - start_time

    # Remove temporary directory
    os.chdir(owd)
    tmp_dir.cleanup()

    # Determine if outlier or not
    if failed is not True:
        covariance = np.array([[6.89334331, 2.95163162], [2.95163162, 5.19778598]])
        location = np.array([-0.44574862, 84.27628921])
        precision = np.array([[0.19167295, -0.10884403], [-0.10884403, 0.25419813]])
        offset = -566.1743358820811
        cov = EllipticEnvelope()
        cov.covariance_ = covariance
        cov.location_ = location
        cov.precision_ = precision
        cov.offset_ = offset        
        dist = cov.mahalanobis(np.array([reaction_energy, activation_energy]).reshape(1, -1))
        if dist > 500:
            smiles_stereo = None
            activation_energy = 10**4 # None
            reaction_energy = 10**4 # None
            sa_score = 10**4 # None

    return smiles_stereo, activation_energy, reaction_energy, sa_score, run_time


def substructure_preserver(mol):
    """
    Check for substructure violates
    Return True: contains a substructure violation
    Return False: No substructure violation
    """        
    
    mol = rdkit.Chem.rdmolops.AddHs(mol) # Note: Hydrogens need to be added for the substructure code to work!
    
    if mol.HasSubstructMatch(rdkit.Chem.MolFromSmarts('[H][C@@]1(*)[C@;R2](*)2[C@@]34[C@;R2]5(*)[C;R1](*)=[C;R1](*)[C@;R2](*)([*;R2]5)[C@@]3([*;R1]4)[C@](*)([*;R2]2)[C@@;R1]1([*])[H]')) == True:
        return True # Has substructure! 
    else: 
        return False # Molecule does not have substructure!
    
    
def substructure_violations(mol):
    """
    Check for substructure violates
    Return True: contains a substructure violation
    Return False: No substructure violation
    """
    violation = False
    forbidden_fragments = ['[C-]', '[S-]', '[O-]', '[N-]', '[*+]', '[*-]' '[PH]', '[pH]', '[N&X5]', '*=[S,s;!R]', '[S&X3]', '[S&X4]', '[S&X5]', '[S&X6]', '[P,p]', '[B,b,N,n,O,o,S,s]~[F,Cl,Br,I]', '*=*=*', '*#*', '[O,o,S,s]~[O,o,S,s]', '[N,n,O,o,S,s]~[N,n,O,o,S,s]~[N,n,O,o,S,s]', '[N,n,O,o,S,s]~[N,n,O,o,S,s]~[C,c]=,:[O,o,S,s,N,n;!R]', '*=N-[*;!R]', '*~[N,n,O,o,S,s]-[N,n,O,o,S,s;!R]']
    for ni in range(len(forbidden_fragments)):
        
        if mol.HasSubstructMatch(Chem.MolFromSmarts(forbidden_fragments[ni])) == True:
            violation = True
            break
        else:
            continue

    return violation

    
def get_properties(smi: str, verbose: bool = False, n_procs: int = 1): 
    '''
    Return fitness functions for the design of OPV molecules.

    Args:
        smile: `str` representing molecule
        verbose: `bool` turn on print statements for debugging
        n_procs: `int` number of processes

    Returns:
        Ea: `float` activation energy with group constraints
        Er: `float` reaction energy with group constraints
        sum_Ea_Er: `float` activation + reaction with group and SAS constraints
        diff_Ea_Er: `float` reaction - activation with group and SAS constraints
    '''

    # Check if substructure is present & there are no fragment violations: 
    mol = Chem.MolFromSmiles(smi)
    call_ = substructure_preserver(mol) # False = =Bad!
    substr_vio = substructure_violations(mol) # True == Bad!
    if substr_vio == True or call_ == False: 
        # return -10**4, -10**4, 10**4
        return 10**4, 10**4, 10**4, 10**4
    
    # Normal calculation can procedd (no fitness violations): 
    smiles_stereo, activation_energy, reaction_energy, sa_score, run_time = run_reaction(smi, n_procs=n_procs)

    # assign values
    Ea = activation_energy
    Er = reaction_energy
    sum_Ea_Er = activation_energy + reaction_energy if sa_score <= 6.0 else 10**4
    diff_Ea_Er = -activation_energy + reaction_energy if sa_score <= 6.0 else 10**4

    return Ea, Er, sum_Ea_Er, diff_Ea_Er

if __name__ == '__main__':        
    smi = 'OC1=CC2=C(C=CC=C2)C(=C1)C1=CC(O)=CC2=C1C=CN=C2'
    Ea, Er, sum_Ea_Er, diff_Ea_Er = get_properties(smi)
    print(f'Activation energy: {Ea}')
    print(f'Reaction energy: {Er}')
    print(f'Activation + Reaction: {sum_Ea_Er}')
    print(f'- Activation + Reaction: {diff_Ea_Er}')
