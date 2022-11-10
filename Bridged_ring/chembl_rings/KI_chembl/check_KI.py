import psycopg2
import psycopg2.extras
import argparse


parser = argparse.ArgumentParser(description="ref_lib")
parser.add_argument(
    "--mol",
    dest='mol',
    )

def get_target_bio(line):
    molrengno = line.split()[0]
    with psycopg2.connect(dbname="chembl_29", user="user", password="user", host="192.168.54.19") as conn:
        with conn.cursor(cursor_factory=psycopg2.extras.DictCursor) as cursor:
            cursor.execute("select assay_id from activities where molregno=%s;",(molrengno,))
            assay_id = cursor.fetchall()[0][0]
            cursor.execute("select tid from Assays where assay_id=%s;", (assay_id,))
            TID = cursor.fetchall()[0][0]
            cursor.execute(
                "select pref_name from target_dictionary where TID=%s;", (TID,)
            )
            try:
                TargetName = cursor.fetchall()[0][0]
                if "kinase" in TargetName:
                    print(f"{line[0:-1]} |{TargetName}")
                else:
                    pass
            except:
                pass
            return True

get_target_bio(parser.parse_args().mol)
