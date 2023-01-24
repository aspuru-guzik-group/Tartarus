import pytest
from tartarus import pce

def test_pce():
    smi = 'c1ccccc1'
    dipm, gap, lumo, combined, pce_pcbm_sas, pce_pcdtbt_sas = pce.get_properties(smi)
    print(f'Dipole moment: {dipm}')
    print(f'HOMO-LUMO Gap: {gap}')
    print(f'LUMO: {lumo}')
    print(f'Combined obj: {combined}')
    print(f'PCE_PCBM - SAS: {pce_pcbm_sas}')
    print(f'PCE_PCDTBT - SAS: {pce_pcdtbt_sas}')
    return True