## The following script extracts some information from blast results. The input is the path of the directory of blast results for each 
## spcies (BOT,CONGO,DAR,ERAWAN,GYROMONAS,Hexamita,IT1,MACHU_PICCU,MIS2C,PIG,SOOS4,TRIMITUS,VLADA7), and the ouput is a pickled data frame with three columns of
## 'qseqid', 'sseqid', 'evalue' for each species. For example the output for BOT is a data frame saved as 'BOT_EUK_df.pkl'.
## The script can also be used to create pickled data frames for each species from the prokaryotic blast results. To create pickled data frames for the prokaryotic 
## blast results make the following changes:
## uncomment the following line:
## #all_blast_files = glob.glob('%s/*_PROK.blastp' %species) #for prokaryotes 
## #final_df.to_pickle("%s_PROK_df.pkl" %os.path.basename(os.path.normpath(species)))
## and comment the following line:
## final_df.to_pickle("%s_EUK_df.pkl" %os.path.basename(os.path.normpath(species)))
## all_blast_files = glob.glob('%s/*_EUK.blastp' %species) #for eukaryotes 

import pandas as pd
import os
import glob
import pickle

#directory contains the path to folders containg the blast results for each species
directory = ['/home/users/fardokht/blast_results_all/blast_results_ALL/BOT', '/home/users/fardokht/blast_results_all/blast_results_ALL/CONGO',
             '/home/users/fardokht/blast_results_all/blast_results_ALL/DAR', '/home/users/fardokht/blast_results_all/blast_results_ALL/ERAWAN',
             '/home/users/fardokht/blast_results_all/blast_results_ALL/GYROMONAS', '/home/users/fardokht/blast_results_all/blast_results_ALL/Hexamita',
             '/home/users/fardokht/blast_results_all/blast_results_ALL/IT1', '/home/users/fardokht/blast_results_all/blast_results_ALL/MACHU_PICCU'
             '/home/users/fardokht/blast_results_all/blast_results_ALL/MIS2C', '/home/users/fardokht/blast_results_all/blast_results_ALL/PIG',
             '/home/users/fardokht/blast_results_all/blast_results_ALL/SOOS4', '/home/users/fardokht/blast_results_all/blast_results_ALL/TRIMITUS',
             '/home/users/fardokht/blast_results_all/blast_results_ALL/VLADA7']


for species in directory:
    all_blast_files = glob.glob('%s/*_EUK.blastp' %species) #for eukaryotes 
    #all_blast_files = glob.glob('%s/*_PROK.blastp' %species) #for prokaryotes 
    full_blast_files = []
    for i in all_blast_files:  # if file is empty, remove the file form the list s
        if (os.stat(i).st_size != 0):
            full_blast_files.append(i)
            #
    qseqid_all = []
    sseqid_all = []
    evalue_all = []
    for i in full_blast_files:
        df = pd.read_table(i, header=None)
        colnames = list(range(1, len(df.columns)+1))
        colnames = [str(i) for i in colnames]
        colnames = ['col' + s for s in colnames]
        #name column
        df.columns = colnames
        qseqid_all.extend(list(df['col1']))
        sseqid_all.extend(list(df['col3']))
        evalue_all.extend(list(df['col12']))
    final_df = pd.DataFrame({'qseqid': qseqid_all,
                             'sseqid': sseqid_all,
                             'evalue': evalue_all})
    final_df.to_pickle("%s_EUK_df.pkl" %os.path.basename(os.path.normpath(species)))
    #final_df.to_pickle("%s_PROK_df.pkl" %os.path.basename(os.path.normpath(species)))
