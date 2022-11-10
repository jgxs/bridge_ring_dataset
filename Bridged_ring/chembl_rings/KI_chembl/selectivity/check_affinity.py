import psycopg2
import psycopg2.extras
# https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/schema_documentation.html
# key value table



all_results = []
with psycopg2.connect(dbname="chembl_29", user="user", password="user", host="192.168.54.19") as conn:
    with conn.cursor(cursor_factory=psycopg2.extras.DictCursor) as cursor:
        molrego = "1293757"
        # cursor.execute("select STANDARD_VALUE from ACTIVITIES where MOLREGNO = %s;",(molrego,))
        # cursor.execute("select * from ACTIVITIES where MOLREGNO = %s;",(molrego,))
        cursor.execute("select ASSAY_ID from ACTIVITIES where MOLREGNO = %s;",(molrego,))
        results = cursor.fetchall()
        out_result = []
        for item in results:
            assay_id = item[0]
            cursor.execute("select TID from ASSAYS where ASSAY_ID = %s;",(assay_id,))
            TID = cursor.fetchall()[0][0]
            cursor.execute(
                "select pref_name from target_dictionary where TID=%s;", (TID,)
            )
            pref_name = cursor.fetchall()[0][0]
            if 'kinase' in pref_name:
                cursor.execute("select STANDARD_VALUE from ACTIVITIES where ASSAY_ID = %s;",(assay_id,))
                bio = cursor.fetchall()[0][0]
                cursor.execute("select STANDARD_TYPE from ACTIVITIES where ASSAY_ID = %s;",(assay_id,))
                assay_type = cursor.fetchall()[0][0]
                print(f"{assay_type:<12} {bio}")
                if assay_type == "IC50":
                    out_result.append((bio,pref_name,assay_type))
                else:
                    pass
        if len(out_result) > 2:
            out_result = sorted(out_result, key = lambda x:x[0] )
            print(out_result)
            print(f"{molrego:>8} {abs(out_result[-1][0]/out_result[0][0]):>.2f} |{out_result[0][1]}| {out_result[-1][1]}")
