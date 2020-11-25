import pandas as pd
import pickle
import os

def append_pkl_df(pkl_df ,data_list, col_list):
    pickle_df = pd.read_pickle(pkl_df) #read in an existing pkl df
    for i in data_list:
        temp_df = pd.read_table(i, header=None)
        temp_df = temp_df.iloc[:,col_list] #subset the columns we want: qseqid, sseqid, evalue
        temp_df.columns = pickle_df.columns
        pickle_df = pickle_df.append(temp_df, ignore_index=True) #appending the new dataframe to the existing pickled df
    pickle_df.to_pickle(pkl_df) #writing out the pickled df. the new okl df will overwrite the previous pkl df


    
