# parallelization:
import pandas as pd
import functools

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
    pool = mp.Pool(processes=50)
    ans = pool.map(functools.partial(remove_duplicate_accession, df=df), list(set(df['qseqid'])))
    return ans

df = pd.read_csv('run_result.csv')
result = main(df)

empty_df = pd.DataFrame()
for i in result:
    empty_df = empty_df.append(i, ignore_index = True)

empty_df.to_csv("removed_duplicate_acc.csv", sep='\t', index=False)
