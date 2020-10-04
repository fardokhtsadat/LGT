import pandas as pd
from multiprocessing import Pool
import numpy as np
from math import nan
import os
import glob
import functools


# parse_fasta_headers() parses the headers from a fasta file and stores it in a file. parse_fasta_headers() gets the directory of a fasta file as input and creates
# a file containg the headers. for example parse_fasta_headers('OG0000000.fa') creates a file called 'headers_OG0000000' with the headers in 'OG0000000.fa'.
def parse_fasta_headers(path_to_dir):
    import pathlib
    import os
    path = pathlib.PurePath(path_to_dir)
    afile = path.name
    os.system("grep -e '>' %s >> headers_%s" %(afile, afile.split('.', 1)[0]))
    os.system("sed -i 's/>//g' headers_%s" %(afile.split('.', 1)[0]))

# modify_string(x) makes sure the headers parsed from fasta files are consistent with the qseqids in pickled data frames. modify_string(x) gets the name of a
# file containg the parsed headers from a fasta file, and makes modifications in the names if needed. In this example, modify_string('headers_OG0000000') gets
# the file 'headers_OG0000000' and modifies the headers in this file if needed.
def modify_string(x):
    with open(x) as f:
        content = f.readlines()
    content = [x.strip() for x in content]
    substring = "_alternative"
    for i in range(len(content)):
        if substring in content[i]:
            content[i] = content[i][:content[i].rfind(substring)]
        else:
            content[i] = content[i]
    return content

# First, search_pkl_df() uses the result from modify_string(x) and groups the headers based on the species they belong to.
# this makes the searching of related qseqids in the pickled data frames faster and more efficient since pickled data frame for the same species are not opened
# and closed constantly. Next, search_pkl_df(x) searches for the qseqids in the pickled data frames and create a csv file.
def get_species(x):
    species_list = ['BOT', 'CONGO', 'DAR', 'GYROMONAS', 'Hexamita', 'IT1', 'MACHU_PICCU', 'MIS2C', 'PIG', 'SOOS4','TRIMITUS', 'VLADA7']
    species = {'BOT': [], 'CONGO': [], 'DAR': [], 'GYROMONAS': [], 'Hexamita': [], 'IT1': [], 'MACHU_PICCU': [],'MIS2C': [], 'PIG': [], 'SOOS4': [], 'TRIMITUS': [], 'VLADA7': []}
    for i in range(len(x)):
        if x[i].split('_')[0] in species_list:
            element = x[i].split('_')[0]
            species.setdefault(element, [])
            species[element].append(x[i])
        species = {k: v for k, v in species.items() if v}
    for j in species:
        species[j] = list(set(species[j]))
    return species


#the following function retrives info from pickle dataframe for an orthogroup
def search_pkl_df(x):
    empty_df = pd.DataFrame()
    for i in x:
        pickle_df = '%s_EUK_df.pkl' % (i)
        df = pd.read_pickle(pickle_df)
        for j in range(len(x[i])):
            empty_df = empty_df.append(df.loc[df['qseqid'] == x[i][j]], ignore_index=True)
    #empty_df.to_csv('run_result.csv')
    return empty_df

# remove_duplicate_accession removes duplicate sseqids from the result obtained from search_pkl_df:
# This functioned is called trough main(df).
def remove_duplicate_accession(x, df):
    EUK_df = pd.DataFrame()
    temp = df.loc[df['qseqid'] == x]
    accession_occurrence = (temp['sseqid'].value_counts()).to_dict()
    for j in accession_occurrence:
        if accession_occurrence[j] != 1:
            minimum = temp.loc[temp['sseqid'] == j]['evalue'].min()
            candidate = temp.loc[(temp['sseqid'] == j) & (temp['evalue'] == minimum)]
            candidate = candidate.sample(n=1)
            EUK_df = EUK_df.append(candidate, ignore_index=True)
        else:
            EUK_df = EUK_df.append(temp.loc[temp['sseqid'] == j], ignore_index=True)
    return EUK_df

def main(df):
    import multiprocessing as mp
    pool = mp.Pool(processes=20) #50 processors are used to speed up the process
    ans = pool.map(functools.partial(remove_duplicate_accession, df=df), list(set(df['qseqid'])))
    return ans


# get_hitproportion_meaneval() calculates the proportion of hits and mean e-value for each accession number.
# get_hitproportion_meaneval() gets a dataframe as input and retruns a dataframe with three columns of 'sseqid', 'occurrence', 'mean_eval'
def get_hitproportion_meaneval(df):
    unique_qseqid = len(set(df['qseqid']))
    # get occurance of each accession number
    count = df['sseqid'].value_counts()
    count = count.to_frame()
    # get the mean value for each accession number
    # mean = df.groupby(['col1']).agg([np.average])
    mean = df.groupby('sseqid').mean()
    joined_df = pd.concat([count, mean], axis=1)  # concat the two dataframes
    joined_df.reset_index(inplace=True)  # reset indicies
    joined_df.columns = ['sseqid', 'occurrence', 'mean_eval']
    joined_df['occurrence'] = joined_df['occurrence'].div(unique_qseqid)
    #joined_df.to_csv("acc_hitprop_meaneval.csv", sep='\t', index=False)
    return joined_df

# attention : edit_acc() is renamed to remove_version_number()
# remove_version_number() gets a list of accession numbers and removes the version number from them.
# the accession numbers in mysql db (EUK_DB) are stored without version numbers.
def remove_version_number(alist):
    edited_acc = []
    for i in alist:
        edited_acc.append(i.split('.', 1)[0])
    return edited_acc

# get_taxid_taxonomy() gets a list of accession numbers and it retrieves the taxon ids for each accession number.
# the output of this function is a dictionary with keys being the accession number and the values being a tuple of taxon id and taxonomy.
def get_taxid_taxonomy(x):
    x = tuple(x)
    import mysql.connector
    mydb = mysql.connector.connect(
        host="localhost",
        user="root",
        password="",
        database="EUK_DB")
    mycursor = mydb.cursor()
    mycursor.execute("SELECT accession, taxid, taxonomy FROM acc_taxid_taxonomy WHERE accession IN %s;" %(x,))
    myresult = mycursor.fetchall()
    print(myresult)
    acc_taxid = {taxid[0]: taxid[1:] for taxid in myresult}
    return acc_taxid

accessions = remove_version_number(list(df1['sseqid']))
chunks = [accessions[i:i + 1000] for i in range(0, len(accessions), 1000)]

with Pool(20) as p:
    res = p.map(get_taxid_taxonomy, chunks)
    Glob = {}  # to join the results from the parallel run
    for i in res:
        Glob.update(i)

#######
def assign_taxid_taxonomy(df, Glob):
    df['taxid'] = np.nan  # add a new col
    #df["taxid"] = pd.to_numeric(df["taxid"])
    df['taxonomy'] = np.nan
    df['sseqid'] = df['sseqid'].apply(lambda x: x.split('.', 1)[0])
    for i in Glob:
        df.loc[df['sseqid'] == i, ['taxid', 'taxonomy']] = int(Glob[i][0]), Glob[i][1]
    return df

# find_name() searches for ser-defined name in the taxonomy, and return a df containing sseqid, taxids, taxonomy, name.

def find_name(list_of_names, Glob, df):
    taxid_taxonomy = {}
    for i in Glob:
        taxid_taxonomy[Glob[i][0]] = Glob[i][1].split(';')
    all_names_df = pd.DataFrame()
    for name in list_of_names:
        temp_list = []
        temp_df = pd.DataFrame()
        for j in taxid_taxonomy:
            if name in taxid_taxonomy[j]:
                temp_list.append(j)
        for taxid in temp_list:
            temp = df.loc[df['taxid'] == taxid]
            temp['name'] = np.nan
            temp.loc[temp['taxid'] == taxid, 'name'] = name
            temp_df = temp_df.append(temp)
        all_names_df = all_names_df.append(temp_df)
    return all_names_df


