import pandas as pd
import numpy as np
from math import nan
import os
import glob
import pathlib
import mysql.connector

# parse_fasta_headers() parses the headers from a fasta file and stores it in a file. parse_fasta_headers() gets the directory of a fasta file as input and creates
# a file containg the headers. for example parse_fasta_headers('OG0000000.fa') creates a file called 'headers_OG0000000' with the headers in 'OG0000000.fa'.
def parse_fasta_headers(path_to_dir):
    path = pathlib.PurePath(path_to_dir)
    afile = path.name
    os.system("grep -e '>' %s >> headers_%s" %(afile, afile.split('.', 1)[0]))
    os.system("sed -i 's/>//g' headers_%s" %(afile.split('.', 1)[0]))

# modify_string(x) makes sure the headers parsed from fasta files are consistent with the qseqids in pickled data frames. modify_string(x) gets the name of a
# file containg the parsed headers from a fasta file, and makes modifications in the names if needed. In this example, modify_string('headers_OG0000000') gets
# the file 'headers_OG0000000' and modifies the headers in this file if needed.
# input is a file containg qseqids separated by '\n'
# output is a list containing modified qseqids
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
    species_list = ['BOT', 'CONGO', 'DAR', 'GYROMONAS', 'Hexamita', 'IT1', 'MACHU_PICCU', 'MIS2C', 'PIG', 'SOOS4','TRIMITUS', 'VLADA7', 'ERAWAN']
    species = {'BOT': [], 'CONGO': [], 'DAR': [], 'GYROMONAS': [], 'Hexamita': [], 'IT1': [], 'MACHU_PICCU': [],'MIS2C': [], 'PIG': [], 'SOOS4': [], 'TRIMITUS': [], 'VLADA7': [], 'ERAWAN': []}
    for i in range(len(x)):
        if x[i].split('_')[0] in species_list:
            element = x[i].split('_')[0]
            species.setdefault(element, [])
            species[element].append(x[i])
        species = {k: v for k, v in species.items() if v}
    for j in species:
        species[j] = list(set(species[j]))
    return species


#the following function retrives info from pickle dataframe for an orthogroup, this function is called from main_search_pkl_df
def search_pkl_df(x):
    alist = []
    for i in x:
        pickle_df = '/home/users/fardokht/EUK_pickled_df/%s_EUK_df.pkl' % (i)
        df = pd.read_pickle(pickle_df) #read in df
        local_df = pd.DataFrame() #an empty df
        local_df['qseqid'] = x[i] #local_df is equal to the values of the key 'i'
        alist.append(local_df.merge(df, left_on='qseqid', right_on='qseqid'))
    merged_df = pd.DataFrame()
    for j in alist:  # the for-loop merges the results retrieved from multiple processors
        merged_df = merged_df.append(j, ignore_index=True)
    merged_df['sseqid'] = merged_df['sseqid'].apply(lambda x: x.split('.', 1)[0])
    return merged_df

# remove_duplicate_accession removes duplicate sseqids from the result obtained from search_pkl_df:
def remove_duplicate_accession(df):
    df = df.groupby(['qseqid','sseqid'], as_index=False)['evalue'].min()
    return df

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
    mydb = mysql.connector.connect(
        host="mole",
        user="fardokht",
        password="21@hgzIyeu*Ahs",
        database="EUK_PROK_DB")
    mycursor = mydb.cursor()
    mycursor.execute("SELECT accession, taxonomy FROM accession_taxid_taxonomy WHERE accession IN %s;" %(x,))
    myresult = mycursor.fetchall()
    acc_taxid = {taxid[0]: taxid[1] for taxid in myresult}
    mycursor.close()
    mydb.close()
    return acc_taxid

#######
def assign_taxid_taxonomy(df, Glob):
    df['taxonomy'] = np.nan
    #df['sseqid'] = df['sseqid'].apply(lambda x: x.split('.', 1)[0])
    df['taxonomy'] = df["sseqid"].map(Glob)
    return df

# find_name() searches for ser-defined name in the taxonomy, and return a df containing sseqid, taxids, taxonomy, name.
def find_name(list_of_names, Glob, df):
    acc_taxonomy = {}
    for i in Glob:
        acc_taxonomy[i] = Glob[i].split(';')
    all_names_df = pd.DataFrame()
    for name in list_of_names:
        print(name)
        temp_list = []
        for j in acc_taxonomy:
            if name in acc_taxonomy[j]:
                temp_list.append(j)
        temp = df[df.sseqid.isin(temp_list)]
        temp['name'] = name
        all_names_df = all_names_df.append(temp)
    return all_names_df

def sort_and_select(df):
    df = df.sort_values(['occurrence', 'mean_eval'], ascending=[False, True])
    # best 3 elements for each taxon id
    top_df = df.groupby('name').head(3) #we get the top 3
    #top_df = top_df.dropna() #remove nan values
    # 2 random elements for each taxon id.
    random = pd.concat([df, top_df]).drop_duplicates(keep=False)
    #random = random.dropna()
    unique_names =list(pd.unique(random['name'])) #get unique rank names  ##attention, nan values are included and have to be removed
    #unique_names = [x for x in unique_names if str(x) != 'nan']
    #unique_names[np.isnan(unique_names)] = 0
    selected_random = pd.DataFrame()
    for i in unique_names:
        if len(random[random['name'] == i]) >= 2: #if there are more than two elements for this taxid
            selected_random = selected_random.append(random[random['name'] == i].sample(n=2, replace=False))
        else:
            selected_random = selected_random.append(random[random['name'] == i].sample(n=1, replace=False))
        acc_removed = list(selected_random['sseqid'])
        random = random[~(random['sseqid'].isin(acc_removed))] #removes accession number which are already sampled. this is because one accession number can be found for different names
        top_df.append(selected_random).to_csv('candidate_accession_numbers.csv')
        
        
def wrapper(afile_qseqids, list_of_names):
    headers = modify_string(afile_qseqids)
    print('input is modified')
    species = get_species(headers)
    df = search_pkl_df(species)
    print('accession numbers and e-values are retrieved')
    #
    df = remove_duplicate_accession(df)
    print('duplicated accession numbers are removed')
    df = get_hitproportion_meaneval(merged_df)
    #
    print('hit proportions and average e-values are caluclated')
    accessions = list(df['sseqid'])
    chunks = [accessions[i:i + 15000] for i in range(0, len(accessions), 15000)]
    #
    res = []
    for i in chunks:
        res.append(get_taxid_taxonomy(i))
    Glob = {}  # to join the results from the parallel run
    for j in res:
        Glob.update(j)  # Glob contains {accession:taxonomy}
    #
    print('taxid and taxonomy is obtained')
    df1 = assign_taxid_taxonomy(df, Glob)
    df2 = find_name(list_of_names, Glob, df1)
    print('searching for names ...')
    sort_and_select(df2)
    print('a csv file with candidate accession numbers is created')
