# the following script fills in a mysql db
# input 1: a file or a directory containing accession number, taxon id, and taxonomy separated by ','.
# input 2: fill_db() also requires the password to access mysql.
# batch size by default is 300.

def fill_db(afile, password):
    import mysql.connector as mysql
    results = []
    with open(afile) as data:
        for line in data.read().split("\n"):
            if len(line.split(",")) < 2:
                continue
            else:
                accession = line.split(",")[0]
                taxon_id = line.split(",")[1]
                taxonomy = line.split(",")[2]
                results.append([accession, taxon_id, taxonomy])

    total = len(results)
    print("Total Number Of Elements is %s." % (total))

    db = mysql.connect(
        host="localhost",
        user="root",
        passwd=password,
        database="EUK_DB"
    )
    cursor = db.cursor()

    batch = []
    counter = 0
    BATCH_SIZE = 300
    cursor_count = 0

    for i in range(len(results)):
        if counter < BATCH_SIZE:
            batch.append(
                [results[i][0], results[i][1], results[i][2]])  # batch => [[ACC1234.0, 'Some full taxa name'], ...]
            counter += 1
        else:
            batch.append([results[i][0], results[i][1], results[i][2]])
            values = (", ").join([f"('{d[0]}', '{d[1]}', '{d[2]}')" for d in batch])
            sql = f"INSERT INTO accession_taxonomy(accession_number, taxon_id, taxonomy) VALUES {values}"
            try:
                cursor.execute(sql)
                db.commit()
            except Exception as exception:
                print(exception)
                print("An error happened in the first try to commit ...")
                try:
                    cursor.execute(sql)
                    db.commit()
                    print("The error is solved ... ")
                except Exception as exception:
                    for i in range(len(batch)):
                        accession = batch[i][0]
                        taxon_id = batch[i][1]
                        file = open("LeftOvers.txt", "a")
                        file.write(accession + " ")
                        file.write(taxon_id + '\n')
                        file.close()
                        print("A file called LeftOvers is created ... ")

            print(cursor.rowcount, "Records Inserted")
            cursor_count += cursor.rowcount
            counter = 0
            batch = []
    else:
        if batch:
            values = (", ").join([f"('{d[0]}', '{d[1]}', '{d[2]}')" for d in batch])
            sql = f"INSERT INTO accession_taxonomy(accession_number, taxon_id, taxonomy) VALUES {values}"

            try:
                cursor.execute(sql)
                db.commit()
            except Exception as exception:
                print(exception)
                print("An error happened in the first try to commit the last batch ...")
                try:
                    cursor.execute(sql)
                    db.commit()
                    print("The error is solved ... ")
                except Exception as exception:
                    for i in range(len(batch)):
                        accession = batch[i][0]
                        taxon_id = batch[i][1]
                        file = open("LeftOvers.txt", "a")
                        file.write(accession + " ")
                        file.write(taxon_id + '\n')
                        file.close()
                        print("A file called LeftOvers is created ... ")
            print(cursor.rowcount, "Records Inserted")
            cursor_count += cursor.rowcount

    print("Total Number Of %s Rows Has Been Added." % (cursor_count))
    db.close()
