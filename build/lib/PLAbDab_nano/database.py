import os
import json
import pandas as pd
import torch

from PLAbDab_nano.sequence_search import SequenceSearch
from PLAbDab_nano.structure_search import StructureSearch
from PLAbDab_nano.database_search.util import add_url_to_data

class PLAbDab_nano(SequenceSearch, StructureSearch):
    """
    The main class for searching data from the nanobody extension of the Patent and Literature Antibody Database (PLAbDab-nano). 
    ...
    Attributes
    ----------
    path_to_db : str
        Path to local PLAbDab-nano.
    n_jobs : int
        Number of threads to use (default 5).

    Methods
    -------
    vhh_seq_search()
        Search PLAbDab-nano with KA-Search for sequences most similar to the query.
    vnar_seq_search()
        Search PLAbDab-nano with BLAST for sequences most similar to the query.
    structure_search()
        Search PLAbDab-nano with NBB2 for structures most similar to the query.
    get_db()
        Retrieve data from PLAbDab-nano.
    column_search()
        Search columns in PLAbDab-nano for specific terms.
    """
    
    def __init__(
        self, 
        path_to_db, 
        n_jobs = 5, 
        **kwargs
    ):
        super().__init__()
        self.config = kwargs
        self.path_to_db = path_to_db
        self.__check_config()
        
        if not os.path.exists(os.path.join(self.path_to_db, "all_sequences.csv.gz")):
            raise ValueError(f"Provided path_to_db ({self.path_to_db}) does not contain an all_sequences.csv.gz file.")
            
        torch.set_num_threads(n_jobs)
        self.all_sequences = pd.read_csv(os.path.join(self.path_to_db, "all_sequences.csv.gz"))
        self.all_sequences = self.all_sequences.fillna({"reference_title": "", "mentioned_antigens": ""})

    def __check_config(self):
        config_file = os.path.join(self.path_to_db, "config.json")

        if os.path.exists(config_file):
            with open(config_file, 'r') as fp:
                config = json.load(fp)

            for key in config:
                if key in self.config:
                    assert config[key] == self.config[key], f"Database was initialised with a different configuration:\n{self.config}"
                else:
                    self.config[key] = config[key]

    def get_db(self, unpaired=True, url=True):

        if unpaired:
            if not url:
                return self.all_sequences
            else:
                return add_url_to_data(self.all_sequences)


    def sequence_plus_structure_search(self, seqs, rmsd_cutoff = 1.25, seq_identity_cutoff = 0.8, filename = "temp_structure.pdb", url=True):
        cdr_seq_search = self.vhh_seq_search(seqs, keep_best_n=100000, regions = ['cdrs'],length_matched=[True], seq_identity_cutoff=seq_identity_cutoff, url=url)
        struc_search = self.structure_search(seqs, rmsd_cutoff = rmsd_cutoff, filename = filename, url = False)
        output = cdr_seq_search[cdr_seq_search.ID.isin(struc_search.ID)]
        output["rmsd"] = struc_search.set_index("ID").loc[output.ID].rmsd.values

        return output.sort_values("rmsd").reset_index(drop=True)


    def column_search(self, term = "antibod", unpaired=True, column=None, url=True, case = False):
        if unpaired:
            df = self.all_sequences
            if column is None:
                column = "reference_title"

        output = df[df[column].str.contains(term, case=case, na=False)].reset_index(drop=True)
        if url and unpaired:
            output = add_url_to_data(output)

        return output
