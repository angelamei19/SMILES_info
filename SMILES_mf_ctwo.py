import shelve
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem as Chem
import multiprocessing
from tqdm import tqdm
import pandas as pd
import pickle
import datamol as dm


class MoleculeAnalyzer:
    def __init__(self, data_file, log_file = 'fp_error_two.txt'):
        self.data_file = data_file
        self.loaded_data = None
        self.log_file = log_file
        self.list_of_all_mols = []
        self.list_of_all_clusters = []
        self.smiles_list = []
    
    
    def load_data(self):
        # Read TSV file into a pandas DataFrame
        self.loaded_data = pd.read_csv(self.data_file, sep='\t')
        self.loaded_data.columns = ['smiles']
        """ with open('sample.pkl', 'rb') as f:
            self.loaded_data = pickle.load(f) """
        self.loaded_data = self.loaded_data.reset_index(drop = True)
        print(self.loaded_data.head(10))
        print(type(self.loaded_data))
        print(len(self.loaded_data))

    
    def analyze_molecules(self):
        # Get list of SMILES strings to process
        self.smiles_list = self.loaded_data['smiles'].tolist()

        # Use multiprocessing.Pool for parallel execution
        with multiprocessing.Pool(25) as pool:
            self.list_of_all_mols = list(tqdm(pool.imap(self.process_smiles, self.smiles_list), total=len(self.smiles_list)))
    
    def process_smiles(self, smiles):
        try:
            # Process each SMILES string independently
            SMILES_mol = Chem.MolFromSmiles(smiles)
            return SMILES_mol
        
        except Exception as e:
            self.log_error(smiles, str(e))
            return None, None
    
    def save_to_shelve(self, shelve_file):
        with shelve.open(shelve_file, 'c') as db:
            for index, row in enumerate(self.list_of_all_clusters):
                db[str(index)] = row

    #Define clustering setup
    def ClusterFps(self):
        with tqdm(total=len(self.list_of_all_mols), desc="Clustering molecules") as pbar:
            clusters, mol_clusters = dm.cluster_mols(self.list_of_all_mols, cutoff=0.2)
            pbar.update(len(self.list_of_all_mols))
        for tup in clusters:
            self.this_cluster = []
            for elem in tup:
                self.this_cluster.append(self.smiles_list[elem])
            self.list_of_all_clusters.append(self.this_cluster)
  

if __name__ == "__main__":
    analyzer = MoleculeAnalyzer('All_SMILES_no_duplicates.tsv')
    #analyzer = MoleculeAnalyzer('sample.pkl')
    analyzer.load_data()
    analyzer.analyze_molecules()
    analyzer.ClusterFps()
    analyzer.save_to_shelve('all_sample_cluster_0.2.shelve')
