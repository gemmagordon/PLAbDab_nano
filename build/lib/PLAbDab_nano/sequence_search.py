import os, time
import numpy as np
import pandas as pd

from kasearch import AlignSequences, SearchDB
from PLAbDab_nano.database_search.util import add_url_to_data
from PLAbDab_nano.database_search.shark_sequence_search import vnar_search, get_metadata


minimum_size_of_results = 1_000_000


class SequenceSearch:
    """
    The main class for searching PLAbDab-nano with KA-Search for sequences most similar to the query. 
    ...
    Methods
    -------
    vhh_seq_search()
        Search PLAbDab-nano with KA-Search for sequences most similar to the VHH query.
    vnar_seq_search()
        Search PLAbDab-nano with BLAST for sequences most similar to the VNAR query.
    """

    def vhh_seq_search(
        self, 
        seqs: dict, 
        keep_best_n: int = 10, 
        seq_identity_cutoff: float = 0.0,
        sort_by: str = 'average',
        regions=['whole'], 
        length_matched=[False],  
        allowed_species=['Human', 'Mouse'], 
        url = False,
        n_jobs=5,
    ):
        """Search database with KA-Search for sequences most similar to the queries.
        
        Parameters
        ----------
        seqs : dict of str
            Heavy chain to search with. If only one chain, search unpaired data.
        keep_best_n : int
            Number of closest matches to return (default is 10).
        seq_identity_cutoff : float
            Lowest sequence identity returned (default is 0.0).
        sort_by : str
            Sort returned sequences by heavy, light or average (default is average).
        regions : list of str
            Region is search across (default is ['whole']).
        length_matched : list of bool
            If search only returns regions with the same length as query (default is [False]).  
        allowed_species : list of str
            Which species to search against (default is ['Human', 'Mouse']).
        url : bool
            If return results with additional column containing the url for the data (default is False)
        n_jobs : int
            Number of threads to use.
            
        """ 

        assert sort_by in ['heavy', 'average'], f"{sort_by} is invalid, use 'heavy' or 'average'"
        
        if len(seqs) == 1: 
            seq = seqs["H"] if "H" in seqs else seqs["L"]
            out = search_unpaired(
                seq, 
                os.path.join(self.path_to_db, "kasearch_db"), 
                regions=regions, 
                length_matched=length_matched,  
                allowed_species=allowed_species, 
                keep_best_n=keep_best_n, 
                seq_identity_cutoff=seq_identity_cutoff,
                n_jobs=n_jobs,
            )
            if url:
                return add_url_to_data(out)
            else:
                return out
        
        else:
            assert False, f"{seqs} needs to be a dict of a single chain."


    def vnar_seq_search(
        self,
        seq,
        data_directory, 
        keep_best_n: int = 10, 
        seq_identity_cutoff: float = 0.0, 
        regions='whole', 
        length_matched=False,
        url=False
    ):
        '''
        Create a BLAST database of VNAR sequences and search against it for most similar sequences.
        '''
        
        if type(seq) == str: 
            out = vnar_search(
                seq, 
                data_directory=data_directory,
                keep_best_n=keep_best_n,
                seq_identity_cutoff=seq_identity_cutoff,
                regions=regions, 
                length_matched=length_matched, 
                url=url 
            )
            return out
        
        else:
            assert False, f"{seq} needs to be a string of a single sequence."


def search_unpaired(
    seq,
    db_path, 
    regions=['whole'], 
    length_matched=[False], 
    allowed_species=['Human', 'Mouse'], 
    keep_best_n=10, 
    seq_identity_cutoff=0.0,
    n_jobs=1,
):
    
    aligned_seqs = AlignSequences(allowed_species=allowed_species, n_jobs=n_jobs)(seq)
    
    chain_db = SearchDB(database_path=db_path, allowed_chain='Mix', regions=regions, length_matched=length_matched)
    chain_db.search(aligned_seqs[:1], keep_best_n=keep_best_n)
    output = chain_db.get_meta(n_query = 0, n_region = 0, n_sequences = 'all', n_jobs = n_jobs)

    return output[output.Identity >= seq_identity_cutoff].copy()

# Search functions for VNARs in shark_sequence_search