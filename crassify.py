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
    '''Takes a folder of protein fasta files as input, maps the genome to protein ID's and gets the genome length.'''
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
                print(nt_ID)
                print('Genome length:', sum(record_lens)*3)
                record_ids = [record.id for record in records]
                proteins[nt_ID] = record_ids
                genome_lens[nt_ID] = genome_len
            genome_protein_map = pd.DataFrame.from_dict(proteins, orient='index')
            genome_lens = pd.DataFrame.from_dict(genome_lens, orient='index')
            genome_lens.columns = ['genome_length']
            print(genome_lens.head())
            print(genome_protein_map.head())
            protein_metadata = pd.merge(genome_protein_map, genome_lens, left_index=True, right_index=True)
            protein_metadata.to_csv('protein_dict.csv')
            return protein_metadata
    else:
        print(f'{path_to_prots} is empty or does not contain fasta files.\nCheck -p path/to/proteins.')


def merge_metadata_maps(maps, path_to_metadata):
    # Merge with metadata to get genome length
    metadata = pd.read_csv(path_to_metadata)
    maps = maps.reset_index()
    maps = maps.rename(columns={'index': 'ID'})
    print(metadata.columns)
    # Make 2 new columns to avoid '.' contained in ncbi_acc / ID merge issues
    maps['ID'] = maps['ID'].str.split('.').str[0]
    metadata['ID'] = metadata['id'].str.split('.').str[0]
    #df = pd.merge(maps, metadata[['id','length']], left_on='ncbi_acc', right_on='id', how='left')
    df = pd.merge(maps, metadata[['ID','length']], left_on='ID', right_on='ID', how='left')
    # Transpose df to have columns ncbi_acc, prot_acc, length
    melt = pd.melt(df, id_vars=['ID','length'])
    melt = melt[melt.value.notna()]
    melt = melt[['ID','length','value']] \
        .sort_values(by='ID') \
        .reset_index(drop = True) \
        .rename(columns={'length':'genome_length'})
    melt.to_csv('metadata_melt.csv', index=False)
    return melt

def merge_metadata_matches(matches, metadata):
    print(f'{datetime.now()}: Merging metadata with matches.tsv protein alignments.')
    pbar = tqdm(total=100)
    matches = pd.read_csv(path_to_matches, sep='\t', header=None, low_memory=False)
    matches.columns = ['qseqid','sseqid', "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
    matches = matches.merge(metadata, left_on=['qseqid'], right_on=['value']) \
        .rename(columns={'genome_length':'qseqid_genome_length','ID':'qseqid_genome'}) \
        .drop('value', axis=1)
    pbar.update(45)
    matches = matches.merge(metadata, left_on=['sseqid'], right_on=['value']) \
        .rename(columns={'genome_length':'sseqid_genome_length','ID':'sseqid_genome'}) \
        .drop('value', axis=1)
    pbar.update(45)
    matches = matches[matches.qseqid_genome != matches.sseqid_genome]
    #matches.to_csv('matches_merged_07_06_23.csv')
    pbar.update(10)
    print(matches.head())
    return matches

def calc_dist(df):
    print(f'{datetime.now()}: Calculating distance between genomes.')
    pbar = tqdm(total=100)
    df_missing = df[(df.qseqid_genome_length.isna())|(df.sseqid_genome_length.isna())]
    #df_missing.to_csv('matched_merged_dists_missing_metadata.csv',index=False)
    dists = df[(df.qseqid_genome_length.notna()) & (df.sseqid_genome_length.notna())]
    #dists.to_csv('matched_merged_dists_metadata_nona.csv',index=False)
    df = dists.sort_values(by=['qseqid','sseqid','sstart','send'])
    pbar.update(50)
    s = df.groupby(['qseqid', 'sseqid'])['send'].shift()
    df['length'] = df['sstart'].lt(s) * df['sstart'].sub(s.fillna(0)) + df['length']
    df['conserved_AA_#'] = (df['length']/100)*df['pident']
    df[['conserved_AA_#','qseqid_genome_length','sseqid_genome_length']] = df[['conserved_AA_#','qseqid_genome_length','sseqid_genome_length']].astype(float)
    df_grp = df.groupby(['qseqid_genome', 'sseqid_genome']).agg({'conserved_AA_#':'sum', 
                                                        'qseqid_genome_length':'first',
                                                        'sseqid_genome_length':'first',
                                                        'length':'sum'}).reset_index()
    df_grp.to_csv('matched_merged_dists.csv',index=False)
    df_grp['distance'] = 1-df_grp['conserved_AA_#']/((df_grp['qseqid_genome_length']+df_grp['sseqid_genome_length'])/2)
    df_dist = df_grp.round(6)
    df_dist.to_csv('matched_merged_dists.csv',index=False)
    pbar.update(50)
    return df_dist 

def to_phylip(df_dist):
    '''Converts distances to PHYLIP format.'''
    print(f'{datetime.now()}: Converting distance matrix to PHYLIP format.')
    pbar = tqdm(total=100)
    matrix = df_dist.pivot_table(columns='qseqid_genome', index='sseqid_genome', values='distance').reset_index()
    matrix = matrix.set_index('sseqid_genome')
    #matrix.to_csv('dist_matrix_phylim_raw.dist', header=True, index=True, sep =' ')
    # Check dimensions of matrix
    cols = matrix.columns.to_list()
    indx = matrix.index.to_list()
    cols_missing = list(set(cols) - set(indx))
    idx_missing = list(set(indx) - set(cols))
    pbar.update(50)
    if len(cols)>len(indx):
        # insert extra column values into index
        for col in cols_missing:
            matrix = matrix.reset_index()
            row_to_append = (matrix.shape[1]-1) * [1]
            row_to_append.insert(0, col)
            matrix.loc[matrix.shape[0]] = row_to_append
            matrix = matrix.set_index('sseqid_genome')
    elif len(indx)>len(cols):
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
    matrix.to_csv('dist_matrix_phylim.dist', header=True, index=True, sep =' ')
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
    # parser.add_argument("-i", dest="diamond_matches", required=True,
    #                     help="path to diamond matches.tsv", metavar="MATCHES")
    parser.add_argument("-p", dest="path_to_proteins", required=True,
                        help="path to proteomes.", metavar="PROTEINS")
    # parser.add_argument("-m", dest="path_to_metadata", required=True,
    #                     help="path to metadata.", metavar="METADATA")
    args = parser.parse_args()

    genome_prot_map = map_genome_to_proteome(args.path_to_proteins)
    # maps = merge_metadata_maps(genome_prot_map, args.path_to_metadata)
    # matches = merge_metadata_matches(args.diamond_matches, maps)
    # dists = calc_dist(matches)
    # phylip = to_phylip(dists)
    # tree = plot_tree('/home/administrator/phd/crass_DB/May_2023_317_proteomes/crassify/dist_matrix_phylim.newick', args.path_to_metadata)



