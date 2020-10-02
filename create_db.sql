CREATE DATABASE EUK_DB;

CREATE TABLE accession_taxid_taxonomy (
id int(22) NOT NULL AUTO_INCREMENT,
accession varchar(64) NOT NULL,
taxid varchar(16) NOT NULL,
taxonomy varchar(800) NOT NULL,
primary key (id),
unique key accession_index (accession) ,
unique key taxonomy (accession));


# IMPORT DATA:
LOAD DATA LOCAL INFILE 'afile'
INTO TABLE acc_taxid_taxonomy
FIELDS TERMINATED BY ','
ENCLOSED BY '"'
LINES TERMINATED BY '\n'
(accession, taxid, taxonomy);

