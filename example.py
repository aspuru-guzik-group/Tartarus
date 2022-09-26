from tartarus import pce
from tartarus import tadf
from tartarus import docking
from tartarus import reactivity
import pdb; pdb.set_trace()

# test and print all objecties
dipm, gap, lumo, combined, pce_pcbm_sas, pce_pcdtbt_sas = pce.get_properties('c1sc(-c2[SiH2]c(cc2)-c2ccc(-c3scc4occc34)c3cscc23)c2Cccc12')
print('******* PCE *******')
print(f'Dipole: {dipm}')
print(f'HOMO-LUMO Gap: {gap}')
print(f'LUMO: {lumo}')
print(f'Combined obj: {combined}')
print(f'PCE1: {pce_pcbm_sas}')
print(f'PCE2: {pce_pcdtbt_sas}')
print()

dipm, gap, lumo, combined = pce.get_surrogate_properties('c1sc(-c2[SiH2]c(cc2)-c2ccc(-c3scc4occc34)c3cscc23)c2Cccc12')
print(f'Dipole: {dipm}')
print(f'HOMO-LUMO Gap: {gap}')
print(f'LUMO: {lumo}')
print(f'Combined obj: {combined}')
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

Ea, Er, sum_Ea_Er, diff_Ea_Er  = reactivity.get_properties('CC=CC(C)=CC=CC=CC1CC2CC1C13CC21C1C=CC3C1', 
    n_procs=1)  # set number of processes
print('******* Reactivity *******')
print(f'Activation energy: {Ea}')
print(f'Reaction energy: {Er}')
print(f'Activation + reactivity: {sum_Ea_Er}')
print(f'Reactivity - activation: {diff_Ea_Er}')
print()
