import pandas as pd
from multiprocessing import Pool
import numpy as np
from math import nan
import os
import glob
import functools

# parse_fasta_headers() parses the headers from a fasta file and stores it in a file. parse_fasta_headers() gets name of a fasta file as input and creates
# a file containg the headers. for example parse_fasta_headers('OG0000000.fa') creates a file called 'headers_OG0000000' with the headers in 'OG0000000.fa'.
def parse_fasta_headers(afile):
    import os
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
    pool = mp.Pool(processes=50) #50 processors are used to speed up the process
    ans = pool.map(functools.partial(remove_duplicate_accession, df=df), list(set(df['qseqid'])))
    return ans

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

#--------------get taxid
def edit_acc(alist):
    edited_acc = []
    for i in alist:
        edited_acc.append(i.split('.', 1)[0])
    return edited_acc

def get_taxid(x):
    acc_taxid = {}
    result = os.popen("blastdbcmd -entry %s -outfmt '%%T' -db /Data/nr_old/nr" % (x)).read() #try batch entry
    try:
        acc_taxid[x] = [int(s) for s in result.split() if s.isdigit()][0]
    except:
        #result = os.popen("efetch -db protein -id %s -format docsum | xtract -pattern DocumentSummary -element TaxId" % (x)).read()
        acc_taxid[x] = ''.join([int(s) for s in result.split() if s.isdigit()])
    return acc_taxid

#----------get taxonomic rank
def get_taxonomic_rank(afile):
    import os
    result = os.popen("epost -db taxonomy -input %s | efetch -format xml | xtract -pattern Taxon -first TaxId -block \"LineageEx/Taxon\" -unless Rank -equals \"no rank\" -sep \'\t\' -element Rank,ScientificName" %(afile)).read()
    result = result.split('\n')
    for i in range(len(result)):
        result[i] = result[i].split('\t')
    taxid_acc = {taxid[0]: taxid[1:] for taxid in result}
    empty = ''
    if empty in taxid_acc:
        del taxid_acc[empty]
    return taxid_acc

#-----------------------RUN-------------------------#
if __name__ == '__main__':
    parse_fasta_headers('OG0000000.fa') #parse_fasta_headers parses the headers from a fasta file. ['OG0000000.fa' contains 555 headers]
    headers = modify_string('headers_OG0000000') #modify_string() gets a file containg the headers separated by '\n'.
    df = search_pkl_df(headers)
    result = main(df) #main(df) removes duplicate sseqids from the df. One qseqid can hit the same sseqid (sequence) multiple times. These duplicated sseqids are removed.

    empty_df = pd.DataFrame()
    for i in result: # the for-loop merges the results retrieved from multiple processors
        empty_df = empty_df.append(i, ignore_index = True)

    df1 = get_hitproportion_meaneval(empty_df)

    with Pool(20) as p:
        res = p.map(get_taxid, edit_acc(list(df1['sseqid'])))
        Glob = {} # to join the results from the parallel run
        for i in res:
            Glob.update(i)

    df1['taxid'] = np.nan #add a new col
    for i in range(len(df1['sseqid'])):
        edited = df1['sseqid'][i].split('.', 1)[0]
        df1.loc[df1['sseqid'] == df1['sseqid'][i], 'taxid'] = Glob[edited]

    unique_taxids = list(set(df1['taxid']))
    with open('all_taxids.txt', 'w') as f:
        for item in unique_taxids:
            f.write("%s\n" % item)

    os.system("split -l 199 all_taxids.txt taxid_part")
    files = list(glob.glob(os.path.join('/home/fardokht/','taxid_part*'))) # a list of file with a specifc name #the dirctory muct be changed
    os.system('rm all_taxids.txt')

    taxid_full_rank = {} # this dict contain taxid as keys and a list of taxonomic division as the values
    for i in files:
        taxid_full_rank.update(get_taxonomic_rank(i))

    os.system('rm taxid_part*')

    taxid_rank = {}
    for i in taxid_full_rank:
        if "phylum" in taxid_full_rank[i]:
            for j in range(len(taxid_full_rank[i])):
                if taxid_full_rank[i][j] == "phylum":
                    taxid_rank[i] = [taxid_full_rank[i][j], taxid_full_rank[i][j+1]]
                continue
        else:
            taxid_rank[i] = [taxid_full_rank[i][0], taxid_full_rank[i][1]]

    # add two more columns 'rank', 'name'
    df1['rank'] = np.nan
    df1['name'] = np.nan

    # taxid_rank is a dict {'taxid': ['rank':'name'], ...}
    #the follwoing, iterates through taxid_rank and for equla values with joined_df['taxid'], it adds rank and name to joined_df['rank'] and joined_df['name']
    for i in taxid_rank:
        df1.loc[df1['taxid'] == int(i), ['rank', 'name']] = [taxid_rank[i][0], taxid_rank[i][1]]

    phylum_only = ['Acidobacteria', 'Aquificae', 'Caldiserica', 'Candidatus Cryosericota', 'Calditrichaeota', 'Chrysiogenetes', 'Coprothermobacterota', 'Deferribacteres', 'Dictyoglomi', 'Elusimicrobia', 'Bacteroidetes', 'Balneolaeota', 'Candidatus Kapabacteria', 'Chlorobi', 'Ignavibacteriae', 'Rhodothermaeota', 'candidate division LCP-89', 'candidate division Zixibacteria', 'Candidatus Cloacimonetes', 'Candidatus Fermentibacteria', 'Candidatus Hydrogenedentes', 'Candidatus Kryptonia', 'Candidatus Latescibacteria', 'Candidatus Marinimicrobia', 'Fibrobacteres', 'Gemmatimonadetes', 'Fusobacteria', 'Krumholzibacteriota', 'Candidatus Tectomicrobia', 'Nitrospinae', 'Nitrospirae', 'Proteobacteria', 'Candidatus Abyssubacteria', 'Candidatus Aureabacteria', 'Candidatus Omnitrophica', 'Chlamydiae', 'Kiritimatiellaeota', 'Lentisphaerae', 'Planctomycetes', 'Verrucomicrobia', 'Spirochaetes', 'Synergistetes', 'Abditibacteriota', 'Actinobacteria', 'Armatimonadetes', 'Candidatus Dormibacteraeota', 'Candidatus Eremiobacteraeota', 'Chloroflexi', 'Candidatus Margulisbacteria', 'Candidatus Melainabacteria', 'Cyanobacteria', 'Deinococcus-Thermus', 'Firmicutes', 'Tenericutes', 'Thermodesulfobacteria', 'Thermotogae']

    #now joined_df is sorted based on occurrence descending and evalue ascending:
    #from the sorted df, we will look for elements of phylum_only in joined_df, and for each element, we get the top 3 acc and 2 random acc numbers.
    joined_df = df1.sort_values(['occurrence', 'mean_eval'], ascending=[False, True])
    top_df = joined_df.groupby('name').head(3) #we get the top 3
    top_df = top_df.dropna() #remove nan values

    random = pd.concat([joined_df, top_df]).drop_duplicates(keep=False)
    unique_names =list(pd.unique(random['name'])) #get unique rank names  ##attention, nan values are included and have to be removed
    unique_names = [x for x in unique_names if str(x) != 'nan']
    #unique_names[np.isnan(unique_names)] = 0

    selected_random = pd.DataFrame()
    for i in unique_names:
        if len(random[random['name'] == i]) >= 2:
            selected_random = selected_random.append(random[random['name'] == i].sample(n=2, replace = False))
        else:
            selected_random = selected_random.append(random[random['name'] == i].sample(n=1, replace=False))

    # join top_df and selected_random dfs.
    top_df = top_df.append(selected_random)

    candidates_df = pd.DataFrame()
    for i in phylum_only:
        candidates_df = candidates_df.append(top_df.loc[top_df['name'].str.lower() == i.lower()])
    candidates_df.to_csv('candidates_df.csv')



