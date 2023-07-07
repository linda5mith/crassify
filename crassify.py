from Bio import SeqIO
import pandas as pd
import numpy as np
import argparse
import subprocess
from io import StringIO
from tqdm import tqdm
import os
import re
from datetime import datetime
pd.set_option('display.max_columns', None)

def map_genome_to_proteome(path_to_prots):
    '''Takes a folder of protein fasta files as input, 
    maps the genome to protein ID's and gets the genome length.'''
    if os.path.isdir(path_to_prots):
        prots = list_files(path_to_prots)
        if prots:
            print(f'{datetime.now()}: Mapping genome ID to protein ID\'s.')
            proteins={}
            genome_lens={}
            for seq in tqdm(prots):
                nt_ID = get_basename(seq)
                records = list(SeqIO.parse(seq,'fasta'))
                record_lens = [len(record.seq) for record in records]
                genome_len = sum(record_lens)*3
                record_ids = [record.id for record in records]
                proteins[nt_ID] = record_ids
                genome_lens[nt_ID] = genome_len
            genome_protein_map = pd.DataFrame.from_dict(proteins, orient='index')
            genome_lens = pd.DataFrame.from_dict(genome_lens, orient='index')
            genome_lens.columns = ['genome_length']
            protein_metadata = pd.merge(genome_protein_map, genome_lens, left_index=True, right_index=True)
            protein_metadata = protein_metadata.reset_index()
            protein_metadata = protein_metadata.rename(columns={'index': 'ID'})
            # Transpose df to have columns ncbi_acc, prot_acc, genome_length
            melt = pd.melt(protein_metadata, id_vars=['ID','genome_length'])
            melt = melt[melt.value.notna()]
            melt = melt[['ID','genome_length','value']] \
                .sort_values(by='ID') \
                .reset_index(drop = True)
            melt.to_csv('metadata_melt.csv', index=False)
            return melt
    else:
        print(f'{path_to_prots} is empty or does not contain fasta files.\nCheck -p path/to/proteins.')

def merge_metadata_matches(matches, input_metadata, db_metadata):
    '''Merges protein metadata with diamond output matches.tsv'''
    print(f'{datetime.now()}: Merging metadata with matches.tsv protein alignments.')
    pbar = tqdm(total=100)
    matches = pd.read_csv(matches, sep='\t', header=None, low_memory=False)
    matches.columns = ['qseqid','sseqid', "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
    matches = matches.merge(input_metadata, left_on=['qseqid'], right_on=['value']) \
        .rename(columns={'genome_length':'qseqid_genome_length','ID':'qseqid_ID'}) \
        .drop('value', axis=1)
    pbar.update(45)
    db_metadata = pd.read_csv(db_metadata, low_memory=False)
    matches = matches.merge(db_metadata, left_on=['sseqid'], right_on=['protein_ncbi_acc']) \
        .rename(columns={'genome_length':'sseqid_genome_length','protein_ncbi_acc':'sseqid_protein_acc',
        'genome_ncbi_acc':'sseqid_genome_acc','virus':'sseqid_virus'})
    pbar.update(45)
    matches = matches[matches.qseqid_ID != matches.sseqid_genome_acc]
    pbar.update(10)
    return matches

def calc_dist(df):
    '''Finds best hits (most conserved AA's shared) for each protein in proteome.
    Finds top 5 highest taxonomic matches based on summation of conserved AA's.
    Returns distance matrix based on formula 1-(ORFs_A)+(ORFs_B)/(len(Genome_A+Genome_B)/2)
    '''
    print(f'{datetime.now()}: Calculating distances between genomes.')
    now = datetime.now()
    dt_string = now.strftime("%d/%m/%Y_%H:%M:%S")
    outpath = f'./crassify_{dt_string}'
    pbar = tqdm(total=100)
    try:
        os.mkdir('./crassify_output')
    except Exception as e:
        print(e)
    # Return best ICTV protein hit for qseqid
    df['conserved_AA_#'] = (df['length']/100)*df['pident']
    df_proteins= df.sort_values(by=['qseqid','sseqid','sstart','send'])
    df_proteins = df_proteins.loc[df_proteins.reset_index().groupby(['qseqid'])['conserved_AA_#'].idxmax()]
    df_proteins= df.sort_values(by=['qseqid','conserved_AA_#','sseqid_virus'])
    df_proteins = df_proteins.reset_index(drop=True)
    df_proteins.to_csv('./crassify_output/protein_hits.csv',index=False)
    pbar.update(50)
    df_genome_hits = df.groupby(['qseqid_ID','sseqid_genome_acc','sseqid_virus'])[['conserved_AA_#','length']].sum().reset_index() \
    .sort_values(by=['conserved_AA_#'], ascending=False)
    df_genome_hits = df_genome_hits.groupby('qseqid_ID').head(5).reset_index(drop=True)
    df_meta = df[['sseqid_genome_acc','Realm','Subrealm', 'Kingdom', 'Subkingdom', 'Phylum', 'Subphylum', 'Class',
       'Subclass', 'Order', 'Suborder', 'Family', 'Subfamily', 'Genus','Subgenus', 'Species', 'Exemplar or additional isolate',
       'Virus name(s)', 'Virus name abbreviation(s)','Virus isolate designation', 'Virus GENBANK accession',
       'Virus REFSEQ accession', 'Genome coverage', 'Genome composition','Host source']]
    df_genome_hits = df_genome_hits.merge(df_meta, left_on=['sseqid_genome_acc'], right_on=['sseqid_genome_acc'], how='left') \
        .drop_duplicates(keep='first').reset_index(drop=True)
    df_genome_hits.to_csv('./crassify_output/genome_hits.csv',index=False)
    # calculate distance 
    # group by qseqid to return best protein hit only for that qseqid
    df.to_csv('./crassify_output/df_so_far.csv', index=False)
    df = df.groupby(['qseqid']).agg({'sseqid':'first','conserved_AA_#':'max',
                              'qseqid_ID':'first',
                              'sseqid_genome_acc':'first',
                              'qseqid_genome_length':'first',
                              'sseqid_genome_length':'first',
                              'sseqid_virus':'first',
                              'length':'first'}).reset_index()
    df = df.groupby(['qseqid_ID','sseqid_genome_acc']).agg({'conserved_AA_#':'sum', 
                                                            'qseqid_genome_length':'first',
                                                            'sseqid_genome_length':'first',
                                                            'sseqid_virus':'first',
                                                            'length':'sum'}).reset_index()
    df['distance'] = 1-(df['conserved_AA_#']*3)/((df['qseqid_genome_length']+df['sseqid_genome_length'])/2)
    df = df[df.qseqid_ID != df.sseqid_genome_acc]
    df['distance'] = df['distance'].abs()
    df = df.sort_values(by=['qseqid_ID','distance'])
    df_dist = df.round(6)
    df_dist.to_csv('./crassify_output/dists.csv',index=False)
    pbar.update(50)
    return df_dist 

def to_phylip(df_dist):
    '''Converts distances to PHYLIP format.'''
    print(f'{datetime.now()}: Converting distance matrix to PHYLIP format.')
    pbar = tqdm(total=100)
    matrix = df_dist.pivot_table(columns='qseqid_ID', index='virus', values='distance').reset_index()
    matrix = matrix.set_index('virus')
    #matrix.to_csv('dist_matrix_phylim_raw.dist', header=True, index=True, sep =' ')
    # Check dimensions of matrix
    cols = dist_piv.columns.to_list()
    indx = dist_piv.index.to_list()
    cols_missing = list(set(indx) - set(cols))
    idx_missing = list(set(cols) - set(indx))
    pbar.update(50)
    if len(cols)>len(indx):
        # insert extra column values into index
        for col in cols_missing:
            matrix = matrix.reset_index()
            row_to_append = (matrix.shape[1]-1) * [1]
            row_to_append.insert(0, col)
            matrix.loc[matrix.shape[0]] = row_to_append
            matrix = matrix.set_index('sseqid_genome')
    if len(indx)>len(cols):
        # insert extra index values as columns
        for idx in idx_missing:
            matrix[idx] = 1
    matrix.to_csv('matrix_what_is_happening.csv',index=True, sep=' ')
    # matrix = matrix.set_index('sseqid_genome')
    #sort dataframe by index and column values
    print(matrix.head())
    sorted_index = sorted(matrix.index)
    matrix = matrix.loc[sorted_index]
    matrix = matrix[sorted_index]
    matrix.columns = [''] * len(matrix.columns)
    matrix.index.names = [len(matrix)]
    # (),: - not permitted by fastme in index
    matrix.index = matrix.index.str.replace('\,|\(|\)|\:','_', regex=True)
    # fastme: taxa names must not be longer than 64 characters
    matrix.index = matrix.index.str[:64]
    np.fill_diagonal(matrix.values, 0)
    matrix = matrix.fillna(1) \
        .astype(float) \
        .round(5) \
        .clip(lower = 0)
    matrix.to_csv('./crassify_output/dist_matrix_phylip.dist', header=True, index=True, sep =' ')
    pbar.update(50)
    return matrix

def to_newick(phylip_matrix, path_to_decenttree):
    print(f'{datetime.now()}: Calling DecentTree to convert PHYLIP to Newick format.')
    # ./decenttree -in /data/san/data1/users/linda/crAss_DB/ictv_proposal_seq_translations/dist_matrix_phylim.dist -t NJ -nt 10 -out dist_matrix_NJ.newick
    p = subprocess.call([path_to_decenttree, "-i", phylip_matrix])
    # check p output

def plot_tree(path_to_newick, path_to_metadata):
    # empress tree-plot -t tree.nwk -fm feature-metadata.tsv -o tree-viz
    # feature metadata needs to be tsv for empress
    metadata = read_metadata(path_to_metadata)
    metadata.to_csv('feature-metadata.tsv', sep='\t', index=False)
    return metadata

def is_valid_arg(parser, arg):
    if not os.path.exists(arg):
        parser.error(f'The file {arg} does not exist!')
    else:
        return arg

def list_files(path):
    path_to_files = []
    for root, dir, files in os.walk(path):
        for f in files: 
            if 'fa' in f or 'fasta':
                path_to_files.append(os.path.join(root, f))
    return path_to_files

def get_basename(path_to_file):
    filepath = os.path.splitext(path_to_file)
    base = filepath[0].split('/')[-1]
    return base


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-m", dest="diamond_matches", required=True,
                        help="path to diamond matches.tsv", metavar="MATCHES")
    parser.add_argument("-p", dest="proteins", required=True,
                        help="path to unclassifed proteomes.", metavar="PROTEINS")
    parser.add_argument("-db", dest="db_metadata", required=True,
                        help="path to database metadata.", metavar="METADATA")
    args = parser.parse_args()

    genome_prot_map = map_genome_to_proteome(args.proteins)
    matches = merge_metadata_matches(args.diamond_matches, genome_prot_map, args.db_metadata)
    dists = calc_dist(matches)
    #phylip = to_phylip(dists)
    # tree = plot_tree('/home/administrator/phd/crass_DB/May_2023_317_proteomes/crassify/dist_matrix_phylim.newick', args.path_to_metadata)



