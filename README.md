# LGT
**LGT** is a python package that suggests some candidate accession numbers for a given orthogroup and a user-defined list of species. **LGT** uses pickled data 
frames created by the user from blast results and contain three columns of 'qseqid', 'sseqid', and 'evalue'. For each query sequence id found in the orthogroup 
file provided by the user, the pickled data frames are searched to retreive the subject sequences hit by the query sequence together with their e-values.
