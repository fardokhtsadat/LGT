CREATE DATABASE EUK_DB;

CREATE TABLE acc_taxid_taxonomy (
    id INT NOT NULL AUTO_INCREMENT,
    accession varchar(64) NOT NULL,
    taxid INT NOT NULL,
    taxonomy varchar(800) NOT NULL,
    PRIMARY KEY (id)
);


# IMPORT DATA:
LOAD DATA LOCAL INFILE 'afile'
INTO TABLE acc_taxid_taxonomy
FIELDS TERMINATED BY ','
ENCLOSED BY '"'
LINES TERMINATED BY '\n'
(accession, taxid, taxonomy);

