#convert each dictionary to tab and \n separated txt file
#input is the name of a pickled dictionary (keys are accession numbers and values are taxonid)
#the output is csv file with three columns of accession number, taxonid, taxonomy

#this function needs taxonkit. please see https://github.com/shenwei356/taxonkit

def get_taxonomy(pickled_dict):
    import pickle
    import os
    import pandas as pd
    import numpy as np
    pickle_in = open(pickled_dict,"rb")
    unpickled_dict = pickle.load(pickle_in)
    with open('acc_taxonid.txt', 'w') as f:
        for item in unpickled_dict:
            f.write("%s\t%s\n" %(item, unpickled_dict[item][0]))
    #
    taxids = []
    for i in unpickled_dict:
        taxids.append(unpickled_dict[i][0])
    #
    unique = list(set(taxids))
    with open('only_taxids.txt', 'w') as f:
        for item in unique:
            f.write("%s\n" %(item))
    #
    os.system('taxonkit lineage only_taxids.txt | tee lineage.txt')
    taxa_dict = {}
    with open("lineage.txt") as data:
        for line in data.read().split("\n"):
            if len(line.split(" ")) < 2:
                continue
            else:
                taxid = line.split("\t")[0]
                taxonomy = line.split("\t")[1]
                taxa_dict[taxid] = taxonomy
    #
    #print("reading in acc_taxonid.txt")
    df = pd.read_csv('acc_taxonid.txt', sep = '\t', header=None)
    df['taxonomy'] = np.nan
    for z in taxa_dict:
        df.loc[df[1] == int(z), 'taxonomy'] = taxa_dict[z]
    #
    df.to_csv('%s_df.txt' %(pickled_dict), index=False)
    os.system('rm acc_taxonid.txt')
    os.system('rm only_taxids.txt')
    os.system('rm lineage.txt')



