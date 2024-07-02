import os
import shutil
import subprocess
import pandas as pd
import ast

from PLAbDab_nano.database_search.util import add_url_to_data
from PLAbDab_nano.database_generate.genbank_scrape import get_vnar_annotation

def get_metadata(results, data_directory, url=False):

    # Get metadata - URL link, more info link, organism, targets mentioned, seq identity, reference title 
    vnar_db = pd.read_csv(os.path.join(data_directory,"vnar_sequences.csv.gz"))

    # Filter all VNAR db by the IDs in results
    output = vnar_db.loc[vnar_db['ID'].isin(results['ID'])]
    
    if url:
        output = add_url_to_data(output)

    output['Seq identity (%)'] = results.Identity

    return output


def vnar_search(
    seq,
    data_directory, 
    keep_best_n: int = 10, 
    seq_identity_cutoff: float = 0.0, 
    regions='whole', 
    length_matched=False,
    url=False
):

    '''Sequence search for VNARs

    Parameters
    ----------
    seqs : str
        Sequence to search with. 
    keep_best_n : int
        Number of closest matches to return (default is 10).
    seq_identity_cutoff : float
        Lowest sequence identity returned (default is 0.0).
    regions : str
        Region to search across (default is 'whole'). Either whole seq ('whole') or CDR3 ('cdr3').
    length_matched : bool
        If True search only returns regions with the same length as query (default is False). 
    '''
    # Create/format BLAST database to search against (db will change with updates)
    vnar_db = pd.read_csv(os.path.join(data_directory,"vnar_sequences.csv.gz"))

    # Account for is matching lengths specified:
    if length_matched == True:
        try:
            filtered_db = vnar_db.loc[vnar_db['sequence'].str.len() == len(seq)] # only use seqs of same length in making db to search against
        except:
            if len(filtered_db) == 0:
                print('No results found for that sequence length, defaulting to all lengths.')
            filtered_db = vnar_db 
    else:
        filtered_db = vnar_db 

    # Write seqs to FASTA file
    vnar_search_path = 'PLAbDab_nano/vnar_search/'
    with open(str(vnar_search_path+'db.fa'), 'w') as db:
        for entry in filtered_db.iterrows():
            db.write(str('>'+entry[1]['ID']+'\n'))
            # Account for if region specified - whole seq or CDR3 only
            if regions == 'whole':
                db.write(str(entry[1]['sequence']+'\n'))
            if regions == 'cdr3':
                cdr3 = ast.literal_eval(entry[1]['cdr_sequences'])['CDRH3']
                db.write(str(cdr3+'\n'))
    
    # Create BLAST db with FASTA file 
    make_db = str('makeblastdb -in '+vnar_search_path+'db.fa -dbtype prot -parse_seqids -out '+vnar_search_path+'vnar_db/vnar_db')
    subprocess.run(make_db, shell=True, stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL)
    
    # Turn query into FASTA file
    if regions == 'whole':
        with open(str(vnar_search_path+'query.fa'), 'w') as query:
            query.write(str('>query'+'\n'))
            query.write(str(seq+'\n'))
    if regions == 'cdr3':
        # Pull out CDR3 from query 
        cdr3 = ast.literal_eval(get_vnar_annotation(seq)[0])['CDRH3']
        with open(str(vnar_search_path+'query.fa'), 'w') as query:
            query.write(str('>query'+'\n'))
            query.write(str(cdr3+'\n'))
    
    # Run query against database
    run_search = str('blastp -query '+vnar_search_path+'query.fa -db '+vnar_search_path+'vnar_db/vnar_db -outfmt 6 -out '+vnar_search_path+'search_results')
    subprocess.run(run_search,shell=True, stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL)
    
    # Process results
    results = pd.read_table(os.path.join(vnar_search_path,'search_results'), sep='\t',header=None).loc[:,1:2]
    results.columns = ['ID', 'Identity']
    # Sort in order of seq identity
    results = results.sort_values(by='Identity', ascending=False)

    # Deal with if number requested greater than number of results
    if keep_best_n > len(results):
        keep_best_n = len(results)

    # Remove temp files
    # os.remove(os.path.join(vnar_search_path,'db.fa'))
    # os.remove(os.path.join(vnar_search_path,'query.fa'))
    # os.remove(os.path.join(vnar_search_path,'search_results'))
    # shutil.rmtree(os.path.join(vnar_search_path,'vnar_db'))

    # Get metadata
    output = get_metadata(results, data_directory, url)
    
    # Return based on input identity cutoff and number results requested 
    return output[output['Seq identity (%)'] >= seq_identity_cutoff][0:keep_best_n].copy()

