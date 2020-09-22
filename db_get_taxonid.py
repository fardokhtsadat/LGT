import pickle
import os

with open('/home/fardokht/DB/part_ai') as f:
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
    print(chunk_count)

pickle_out = open("part_ai.pickle", "wb")
pickle.dump(all_acc_taxid, pickle_out)
pickle_out.close()
