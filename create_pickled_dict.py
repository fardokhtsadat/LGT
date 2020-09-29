import os
import pickle
import pandas as pd
import glob

def extract_blast(x):
    df = pd.read_pickle(x)
    temp_dict_accession = {}
    temp_dict_evalue = {}
    # #get unique values
    unique_queryIds = df['qseqid'].unique().tolist()
    for i in range(len(unique_queryIds)):
        evalue = df.loc[df['qseqid'] == unique_queryIds[i]]['evalue'].tolist() # col1 refers to queryIds and col12 refers to the evalues
        accessions = df.loc[df['qseqid'] == unique_queryIds[i]]['sseqid'].tolist() # col 3 refers to seqId = accession_numbers
        temp_dict_accession[unique_queryIds[i]] = accessions
        temp_dict_evalue[unique_queryIds[i]] = evalue
    return temp_dict_accession, temp_dict_evalue


all_pickled_dfs = glob.glob('~/EUK_pickled_df/*_EUK_df.pkl')

for i in all_pickled_dfs:
    accession, evalue = extract_blast(i)
    name = i.split('_')[0]
    pickle_out = open("%s_EUK_accession.pickle" %(name), "wb") # dump and save the dictionary:
    pickle.dump(accession, pickle_out)
    pickle_out.close()
    pickle_out = open("%s_EUK_evalue.pickle" % (name), "wb") # dump and save the dictionary:
    pickle.dump(evalue, pickle_out)
    pickle_out.close()
