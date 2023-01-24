import argparse
import csv
from os import listdir
import multiprocessing as mp
from itertools import repeat
from dataclasses import dataclass, fields

import pandas as pd

from tartarus import pce
from tartarus import tadf
from tartarus import docking
from tartarus import reactivity

@dataclass 
class ResultBase():
    smile: str

@dataclass
class PCEResult(ResultBase):
    dipm: float
    gap: float
    lumo: float
    combined: float
    pce_pcbm_sas: float
    pce_pcdtbt_sas: float

@dataclass
class TADFResult(ResultBase):
    st: float
    osc: float
    combined: float

@dataclass
class DockingResult(ResultBase):
    score_1syh: float
    score_6y2f: float
    score_4lde: float

@dataclass
class ReactivityResult(ResultBase):
    Ea: float
    Er: float
    sum_Ea_Er: float
    diff_Ea_Er: float

class BenchmarkResults():
    def __init__(self, mode, results):
        self.mode = mode 
        self.results = results

    def save(self):
        with open(f'/data/{self.mode}-output.csv', 'w', newline='') as output_file:
            wr = csv.writer(output_file)
            for idx, result in enumerate(self.results):
                if idx == 0:
                    wr.writerow([field.name for field in fields(result)])
                wr.writerow([getattr(result, field.name) for field in fields(result)])

def benchmark_smile(smile, mode):
    if mode == 'pce':
        result = PCEResult(smile, *pce.get_properties(smile))
    elif mode == 'tadf':
        result = TADFResult(smile, *tadf.get_properties(smile))
    elif mode == 'docking':
        result = DockingResult(smile, docking.get_1syh_score(smile), docking.get_6y2f_score(smile), docking.get_4lde_score(smile))
    elif mode == 'reactivity':
        result = ReactivityResult(smile, *reactivity.get_properties(smile))
    else:
        raise ValueError('Invalid mode')

    return result

def benchmark_smiles(smiles, mode, parallel=True):

    if parallel:
        n_procs = mp.cpu_count()
        with mp.Pool(n_procs) as p:
            results = p.starmap(benchmark_smile, zip(smiles, repeat(mode)))
    else:
        results = [benchmark_smile(smile, mode) for smile in smiles]

    return BenchmarkResults(mode, results)

def find_files(path, extension):
    return [f for f in listdir(path) if f.endswith(extension)]

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--data_path', type=str, default="/data", help='Path to data directory')
    parser.add_argument('--output_filename', type=str, default="output.csv", help='Output filename')
    parser.add_argument('--mode', type=str, required=True, help='Mode to run benchmark in')

    args = parser.parse_args()
    return args

if __name__ == "__main__":
    args = get_args()

    for file in find_files('/data', '.csv'):
        df = pd.read_csv(f'/data/{file}')
        smiles = df['smiles'].tolist()
        results = benchmark_smiles(smiles, args.mode)
        results.save()





