## The following script extracts some information from blast results. The input is the path of the directory of blast results for each 
## the ouput is a pickled data frame with three columns of 'qseqid', 'sseqid', 'evalue' for each species.

import pandas as pd, pickle, glob, os, multiprocessing
from multiprocessing import Pool
import functools

def extract_data(directory, cols):
    all_files = []
    for file in glob.glob(directory):
        all_files.append(file)
    full_blast_files = []
    for j in all_files:
        if (os.stat(j).st_size != 0): # check if file is empty
            full_blast_files.append(j)
    qseqid_all = []
    sseqid_all = []
    evalue_all = []
    for z in full_blast_files:
        df = pd.read_table(z, header=None)
        # select columns
        qseqid_all.extend(list(df.iloc[:,cols[0]]))
        sseqid_all.extend(list(df.iloc[:,cols[1]]))
        evalue_all.extend(list(df.iloc[:,cols[2]]))
    final_df = pd.DataFrame({'qseqid': qseqid_all,
                             'sseqid': sseqid_all,
                             'evalue': evalue_all})
    print(multiprocessing.current_process()) #prints the worker thread's number
    return final_df


def extract_data_main(num_cpu, list_of_dir, output_file_name, columns):
    with Pool(num_cpu) as p:
        res = p.map(functools.partial(extract_data, cols = columns), list_of_dir)
    merged_df = pd.DataFrame()
    for i in res:
        merged_df = pd.concat([i, merged_df], ignore_index=True)
    merged_df.to_pickle(output_file_name)
    head, tail = os.path.split(output_file_name)
    print('The pickled dataframe is created and can be found in %s.' % head)
    print('The total number of rows in your pickled dataframe is %i.' % len(merged_df.index))


