import psycopg2
import psycopg2.extras

from rdkit import Chem
patt = Chem.MolFromSmarts("[x3&!a,x4]")
def check_bridgedrings(smi,patt):
    m = Chem.MolFromSmiles(smi)
    if m.HasSubstructMatch(patt):
        return True
    else:
        return False

with psycopg2.connect(dbname="chembl_29", user="user", password="user", host="192.168.54.19") as conn:
    with conn.cursor(cursor_factory=psycopg2.extras.DictCursor) as cursor:
        cursor.execute("select molregno, canonical_smiles from compound_structures;")
        Molrengno = cursor.fetchall()
        with open("bridged_rings_smis.csv","w") as smis:
            for item in Molrengno:
                if check_bridgedrings(item[1],patt):
                   smis.write(f"{item[0]:>8} {item[1]}\n")

# patt = Chem.MolFromSmarts("[x3&!a,x4]")