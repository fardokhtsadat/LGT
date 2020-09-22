# get_taxonid(afile)
# input is a file that contains accession numbers separated by break lines
# output is a pickled dictionary with keys being the accession numbers and the values being taxonids: {accession : ['taxonid']}

def get_taxonid(afile):
    import pickle
    import os
    with open(afile) as f: #afile is 'part_ai'
        accession = f.read().splitlines()
    chunks = [accession[x:x+199] for x in range(0, len(accession), 199)]
    all_acc_taxid = {}
    chunk_count = 0
    for sublist in chunks:
        chunk_count += 1
        with open('taxids_chunk.txt', 'w') as f:
            for item in sublist:
                f.write("%s\n" % item)
        result = os.popen("cat taxids_chunk.txt | epost -db protein | esummary -db protein | xtract -pattern DocumentSummary -element Caption,TaxId").read()
        result = result.split('\n')
        for i in result:
            if i == '':
                result.remove(i)
        for i in range(len(result)):
            result[i] = result[i].split('\t')
        acc_taxid = {taxid[0]: taxid[1:] for taxid in result}
        all_acc_taxid.update(acc_taxid)
        os.system('rm taxids_chunk.txt')
        #print(chunk_count)
    pickle_out = open("%s.pickle" %(afile), "wb")
    pickle.dump(all_acc_taxid, pickle_out)
    pickle_out.close()
