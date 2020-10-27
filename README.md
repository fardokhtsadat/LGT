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

## Installation
The ZIP-File of the package can be downloaded via "Clone or download" as well as the command line:
```markdown
git clone https://github.com/fardokhtsadat/LGT.git
```

## Usage
You can import **LGT** using:
```python
from LGT import *
```
After importing **LGT**, you can call the *wrapper()* function. *wrapper()* receives two arguments; the first argument is the orthogroup file, and the second argument is the list of species. The output of **LGT** is a csv file containg the suggested condidate accession numbers for each species.

## Attention
Please set your working directory to the directory of LGT. You can use the following command to set the working directory:
```python
import os
os.chdir("/path/to/LGT")
```

