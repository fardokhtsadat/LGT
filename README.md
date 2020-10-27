# LGT
**LGT** is a python package that suggests candidate accession numbers for a given orthogroup and a user-defined list of species. **LGT** uses pickled data 
frames created by the user from blast results. These data frames contain three columns of 'qseqid', 'sseqid', and 'evalue'. For each query sequence id found in the orthogroup file provided by the user, the pickled data frames are searched to retreive the subject sequences hit by the query sequence together with their e-values.
Next, the proportion of hits among all query sequence ids and the average e-values are calculated for each subject sequence id. For each subject sequence id, the taxonomics ranks are obtained and searched for the user-defined list of species. Finally, for each species, the three best and two random subject sequence ids are suggested to the user for the creation of phylogenetic trees.

### Author
Please contact me in case of questions, comments, bug reports, etc...

    Author: Fardokhtsadat Mohammadi
    E-Mail: fardokht.fm@gmail.com
    
## Dependencies & System Requirements
**LGT** makes use of several python packages. Here we provide a list of used python packages:
* pandas
* numpy 
* math 
* glob
* mysql.connector
* os