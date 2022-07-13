from tartarus import pce
from tartarus import tadf
from tartarus import docking
from tartarus import reactivity

# test and print all objecties
dipole, hl_gap, lumo, obj, pce_1, pce_2, sas = pce.get_properties('c1sc(-c2[SiH2]c(cc2)-c2ccc(-c3scc4occc34)c3cscc23)c2Cccc12')
print('******* PCE *******')
print(f'Dipole: {dipole}')
print(f'HOMO-LUMO Gap: {hl_gap}')
print(f'LUMO: {lumo}')
print(f'Combined obj: {obj}')
print(f'PCE1: {pce_1}')
print(f'PCE2: {pce_2}')
print(f'SAS: {sas}')
print()

dipole, hl_gap, lumo, obj = pce.get_surrogate_properties('c1sc(-c2[SiH2]c(cc2)-c2ccc(-c3scc4occc34)c3cscc23)c2Cccc12')
print(f'Dipole: {dipole}')
print(f'HOMO-LUMO Gap: {hl_gap}')
print(f'LUMO: {lumo}')
print(f'Combined obj: {obj}')
print()

st, osc, combined = tadf.get_properties('O=C1NC2=C(O1)C1=C(N=CO1)C2=O')
print('******* TADF *******')
print(f'Singlet-triplet: {st}')
print(f'Oscillator strength: {osc}')
print(f'Combined obj: {combined}')
print()

score = docking.get_1syh_score('Cl/C(=C/Cl)/CC(N)C(=O)O')
print('******* Docking *******')
print(f'Docking score: {score}')
print()

activation_energy, reaction_energy, sa_score = reactivity.get_properties('CC=CC(C)=CC=CC=CC1CC2CC1C13CC21C1C=CC3C1', 
    n_procs=1)  # set number of processes
print('******* Reactivity *******')
print(f'Activation energy: {activation_energy}')
print(f'Reaction energy: {reaction_energy}')
print(f'SAS: {sa_score}')
print()
