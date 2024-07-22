import os
import shutil
import subprocess
import pandas as pd
import ast

from PLAbDab_nano.util import add_url_to_data
#from PLAbDab_nano.database_generate.genbank_scrape import get_vnar_annotation
from PLAbDab_nano.region_definitions import *
from PLAbDab_nano.shark_number import number_shark


def get_vnar_annotation(seq):

    '''Get CDR1 & CDR3 sequences for VNARs with language model'''

    try:
        seq = [seq, ''] # seems to bug out if directly put in str seq
        numb = number_shark.number(seq[0:1])[0]

        cdr1 = ''.join([n[1] for n in numb if int(n[0][0]) in reg_def['CDRH1']])#.replace('-','')
        # cdr2 doesn't exist for VNAR, HV2, HV4 instead
        cdr3 = ''.join([n[1] for n in numb if int(n[0][0]) in reg_def['CDRH3']]).replace('<EOS>','') # sometimes last token included in CDR3, these will be filtered out anyway, too short
        cdr_seqs = {'CDRH1':cdr1,'CDRH2':'','CDRH3':cdr3}
        cdr_lens = str(len(cdr1)) + '_|_' + str(len(cdr3))
    
    except Exception as e:
        cdr_seqs, cdr_lens = {},''

    return str(cdr_seqs), cdr_lens


def get_metadata(results, data_directory, url=False):

    # Get metadata - URL link, more info link, organism, targets mentioned, seq identity, reference title 
    vnar_db = pd.read_csv(os.path.join(data_directory,"vnar_sequences.csv.gz"))

    # Filter all VNAR db by the IDs in results
    output = vnar_db.loc[vnar_db['ID'].isin(results['ID'])]
    
    if url:
        output = add_url_to_data(output)

    output['Identity'] = results.Identity

    return output


def vnar_search(
    seq,
    data_directory, 
    vnar_search_path = 'tmp',
    keep_best_n: int = 10, 
    seq_identity_cutoff: float = 0.0, 
    regions='whole', 
    length_matched=[False],
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
    length_matched : list of bool
        If True search only returns regions with the same length as query (default is [False]). 
    url : bool 
        If True return results with additional column containing the url for the data (default is False)
    vnar_search_path : str
        Path to temp directory to write BLAST database and results
    '''
    # Create/format BLAST database to search against (db will change with updates)
    vnar_db = pd.read_csv(os.path.join(data_directory,"vnar_sequences.csv.gz"))

    # Account for is matching lengths specified:
    if length_matched == [True]:
        try:
            filtered_db = vnar_db.loc[vnar_db['sequence'].str.len() == len(seq)] # only use seqs of same length in making db to search against
        except:
            if len(filtered_db) == 0:
                print('No results found for that sequence length, defaulting to all lengths.')
            filtered_db = vnar_db 
    else:
        filtered_db = vnar_db 

    # Write seqs to FASTA file
    # vnar_search_path = 'PLAbDab_nano/PLAbDab_nano/vnar_search/'
    
    # Don't want to delete the directory if it already exists!
    remove_tempdir = True
    if os.path.exists(vnar_search_path):
        remove_tempdir = False
    else:
        os.makedirs(vnar_search_path)
        
    with open(os.path.join(vnar_search_path, 'db.fa'), 'w') as db:
        for entry in filtered_db.iterrows():
            db.write(str('>'+entry[1]['ID']+'\n'))
            # Account for if region specified - whole seq or CDR3 only
            if regions == 'whole':
                db.write(str(entry[1]['sequence']+'\n'))
            if regions == 'cdr3':
                cdr3 = ast.literal_eval(entry[1]['cdr_sequences'])['CDRH3']
                db.write(str(cdr3+'\n'))
    
    # Create BLAST db with FASTA file 
    cmd = [
        'makeblastdb',
        '-in', os.path.join(vnar_search_path, 'db.fa'),
        '-dbtype', 'prot',
        '-parse_seqids',
        '-out', os.path.join(vnar_search_path, 'vnar_db', 'vnar_db')
    ]
    #make_db = str('makeblastdb -in '+vnar_search_path+'db.fa' + '-dbtype prot -parse_seqids -out '+vnar_search_path+'vnar_db/vnar_db')
    p = subprocess.run(cmd, capture_output=True)
    if p.returncode:
        print(p.stderr)
        exit()
        
    # Turn query into FASTA file
    if regions == 'whole':
        with open(os.path.join(vnar_search_path,'query.fa'), 'w') as query:
            query.write(str('>query'+'\n'))
            query.write(str(seq+'\n'))
    if regions == 'cdr3':
        # Pull out CDR3 from query 
        cdr3 = ast.literal_eval(get_vnar_annotation(seq)[0])['CDRH3']
        with open(os.path.join(vnar_search_path,'query.fa'), 'w') as query:
            query.write(str('>query'+'\n'))
            query.write(str(cdr3+'\n'))
    
    # Run query against database
    cmd = [
        'blastp',
        '-query', os.path.join(vnar_search_path, 'query.fa'),
        '-db', os.path.join(vnar_search_path, 'vnar_db', 'vnar_db'),
        '-outfmt', '6',
        '-out', os.path.join(vnar_search_path, 'search_results')
    ]
    #run_search = str('blastp -query '+vnar_search_path+'query.fa -db '+vnar_search_path+'vnar_db/vnar_db -outfmt 6 -out '+vnar_search_path+'search_results')
    p = subprocess.run(cmd, capture_output=True)
    if p.returncode:
        print(p.stderr)
        exit()
            
    # Process results
    results = pd.read_table(os.path.join(vnar_search_path,'search_results'), sep='\t',header=None).loc[:,1:2]
    results.columns = ['ID', 'Identity']
    # Normalise identity between 0-1 to match VHH seq search / for web-app 
    results['Identity'] = results['Identity'].map(lambda Identity: Identity / 100)
    # Sort in order of seq identity
    # results = results.sort_values(by='Identity', ascending=False)

    # Deal with if number requested greater than number of results
    if keep_best_n > len(results):
        keep_best_n = len(results)

    # Remove temp files
    os.remove(os.path.join(vnar_search_path,'db.fa'))
    os.remove(os.path.join(vnar_search_path,'query.fa'))
    os.remove(os.path.join(vnar_search_path,'search_results'))
    
    if not shutil.rmtree.avoids_symlink_attacks: 
        raise RuntimeError("rmtree is not safe on this Platform")        
    shutil.rmtree(os.path.join(vnar_search_path,'vnar_db'))        

    # Get metadata
    output = get_metadata(results, data_directory, url)
    
    # Return based on input identity cutoff and number results requested 
    return output[output['Identity'] >= seq_identity_cutoff].copy().sort_values(by='Identity', ascending=False)[0:keep_best_n]

