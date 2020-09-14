import pandas as pd
import os

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
def search_pkl_df(content):
    import pandas as pd
    species_list = ['BOT', 'CONGO', 'DAR', 'GYROMONAS', 'Hexamita', 'IT1', 'MACHU_PICCU', 'MIS2C', 'PIG', 'SOOS4','TRIMITUS', 'VLADA7']
    species = {'BOT': [], 'CONGO': [], 'DAR': [], 'GYROMONAS': [], 'Hexamita': [], 'IT1': [], 'MACHU_PICCU': [], 'MIS2C': [], 'PIG': [], 'SOOS4': [], 'TRIMITUS': [], 'VLADA7': []}
    for i in range(len(content)):
        if content[i].split('_')[0] in species_list:
            element = content[i].split('_')[0]
            species.setdefault(element, [])
            species[element].append(content[i])
        species = {k: v for k, v in species.items() if v}
    for j in species:
        species[j] = list(set(species[j]))
    empty_df = pd.DataFrame()
    for z in species:
        pickle_df = '%s_EUK_df.pkl' % (z)
        df = pd.read_pickle(pickle_df)
        for h in range(len(species[z])):
            empty_df = empty_df.append(df.loc[df['qseqid'] == species[z][h]], ignore_index=True)
    empty_df.to_csv('run_result.csv')
    
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
    import functools
    pool = mp.Pool(processes=50)
    ans = pool.map(functools.partial(remove_duplicate_accession, df=df), list(set(df['qseqid'])))
    return ans

#--
# edit_acc() removes the version number from an accession number. 
def edit_acc(alist):
    edited_acc = []
    for i in alist:
        edited_acc.append(i.split('.', 1)[0])
    return edited_acc

# get_taxid() gets the taxon id for each accession number.
def get_taxid(x):
    acc_taxid = {}
    result = os.popen("blastdbcmd -entry %s -outfmt '%%T' -db /Data/nr_old/nr" % (x)).read() #try batch entry
    try:
        acc_taxid[x] = [int(s) for s in result.split() if s.isdigit()][0]
    except:
        #result = os.popen("efetch -db protein -id %s -format docsum | xtract -pattern DocumentSummary -element TaxId" % (x)).read()
        acc_taxid[x] = ''.join([int(s) for s in result.split() if s.isdigit()])
    return acc_taxid


#----



if __name__ == '__main__':
    #parse_fasta_headers(afile)
    content = modify_string('headers_OG0000000')
    search_pkl_df(content)
    #
    df = pd.read_csv('run_result.csv')
    result = main(df)
    empty_df = pd.DataFrame()
    for i in result:
        empty_df = empty_df.append(i, ignore_index = True)
    empty_df.to_csv("removed_duplicate_acc.csv", sep='\t', index=False)
    #
    with Pool(20) as p:
        result = p.map(get_taxid, edit_acc(list(joined_df['sseqid'])))
        Glob = {} # to join the results from the parallel run
        for i in result:
            Glob.update(i)
     

    
