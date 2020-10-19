#---required packages---
import pandas as pd
from multiprocessing import Pool
import numpy as np
from math import nan
import os
import glob
import functools
import pathlib
import mysql.connector
#---LGT.py---
import modify_string
import get_species
import search_pkl_df
import remove_duplicate_accession
import get_hitproportion_meaneval
import get_taxid_taxonomy
import assign_taxid_taxonomy
import find_name
import sort_and_select

def wrapper(afile_qseqids, list_of_names, number_of_cpus):
    headers = modify_string(afile_qseqids)
    species = get_species(headers)
    df = search_pkl_df(species)
    #
    with Pool(number_of_cpus) as p:
        all_df = []
        for qseqid, df_qseqid in df.groupby('qseqid'):
            all_df.append(df_qseqid)
        res = p.map(remove_duplicate_accession, all_df)
        merged_df = pd.DataFrame()
        for i in res:  # the for-loop merges the results retrieved from multiple processors
            merged_df = merged_df.append(i, ignore_index=True)
    #
    df = get_hitproportion_meaneval(merged_df)
    #
    accessions = list(df['sseqid'])
    chunks = [accessions[i:i + 10000] for i in range(0, len(accessions), 10000)]
    #
    with Pool(number_of_cpus) as p:
        res = p.map(get_taxid_taxonomy, chunks)
        Glob = {}  # to join the results from the parallel run
        for i in res:
            Glob.update(i)  # Glob contains {accession:taxonomy}
    #
    df1 = assign_taxid_taxonomy(df, Glob)
    df2 = find_name(list_of_names, Glob, df1)
    sort_and_select(df2)
    print('a csv file with candidate accession numbers is created')
    
    

