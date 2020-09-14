import pandas as pd
import os

# parse_fasta_headers() parses the headers from a fasta file and stores it in a file. parse_fasta_headers() gets name of a fasta file as input and creates 
# a file containg the headers. for example parse_fasta_headers('OG0000000.fa') creates a file called 'headers_OG0000000' with the headers in 'OG0000000.fa'.

def parse_fasta_headers(afile):
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


# following groups queries into related species
def get_species(content):
    species_list = ['BOT', 'CONGO', 'DAR', 'GYROMONAS', 'Hexamita', 'IT1', 'MACHU_PICCU', 'MIS2C', 'PIG', 'SOOS4',
                    'TRIMITUS', 'VLADA7']

    species = {'BOT': [], 'CONGO': [], 'DAR': [], 'GYROMONAS': [], 'Hexamita': [], 'IT1': [], 'MACHU_PICCU': [],
               'MIS2C': [], 'PIG': [], 'SOOS4': [], 'TRIMITUS': [], 'VLADA7': []}

    for i in range(len(content)):
        if content[i].split('_')[0] in species_list:
            element = content[i].split('_')[0]
            species.setdefault(element, [])
            species[element].append(content[i])

        species = {k: v for k, v in species.items() if v}

    return species


#the following function retrives info from pickle dataframe for an orthogroup
def hit_acc(x):
    empty_df = pd.DataFrame()
    for i in x:
        pickle_df = '%s_EUK_df.pkl' % (i)
        df = pd.read_pickle(pickle_df)
        for j in range(len(x[i])):
            empty_df = empty_df.append(df.loc[df['qseqid'] == species[i][j]], ignore_index=True)
    empty_df.to_csv('run_result.csv')

content = modify_string('header_OG0000000')
species = get_species(content)
hit_acc(species)
