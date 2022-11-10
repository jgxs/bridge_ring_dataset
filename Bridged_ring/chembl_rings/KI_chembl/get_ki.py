import psycopg2
import psycopg2.extras
# https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/schema_documentation.html
# key value table

from rdkit import Chem
from rdkit.Chem import Descriptors

patt = Chem.MolFromSmarts("[x3]")
patt_spiro = Chem.MolFromSmarts("[x4]")
patt_p_bridge = Chem.MolFromSmarts("[x3&!a]")

with psycopg2.connect(dbname="chembl_29", user="user", password="user", host="192.168.54.19") as conn:
    with conn.cursor(cursor_factory=psycopg2.extras.DictCursor) as cursor:
        with open("/home/chengyj/kinase_work/Bridged_ring/chembl_rings/KI_chembl/ki_chembl.smi") as kis:
            for info in kis:
                molrego = info.split()[0]
                print(molrego)
                try:
                    cursor.execute("select canonical_smiles from compound_structures where molregno=%s;",(molrego,))
                    smis = cursor.fetchall()[0][0]
                    cursor.execute("select assay_id from activities where molregno=%s;",(molrego,))
                    assay_id = cursor.fetchall()[0][0]
                    cursor.execute("select tid from Assays where assay_id=%s;", (assay_id,))
                    TID = cursor.fetchall()[0][0]
                    cursor.execute("select pref_name from target_dictionary where TID=%s;", (TID,))
                    target_name = cursor.fetchall()[0][0]
                    mol = Chem.MolFromSmiles(smis)
                    line = f"{molrego} {smis} | {target_name}\n"
                    with open("data_refined/chembl_ki.csv","a") as ki:
                        ki.write(line)
                    ri = mol.GetRingInfo()
                    largest_ring_size = max((len(r) for r in ri.AtomRings()), default=0)
                    if largest_ring_size > 12 and Descriptors.MolWt(mol) <= 1200:
                        with open("data_refined/bigrings_ki.csv","a") as bigring:
                            bigring.write(line)
                    else: 
                        if mol.HasSubstructMatch(patt_spiro):
                            with open("data_refined/spiro_ki.csv","a") as spiro:
                                spiro.write(line)
                        if mol.HasSubstructMatch(patt_p_bridge):
                            with open("data_refined/bridged_ki.csv", "a") as bridge:
                                bridge.write(line)         
                        if mol.HasSubstructMatch(patt_spiro) and mol.HasSubstructMatch(patt_p_bridge):
                            with open("data_refined/spiro_bridged_ki.csv","a") as sbk:
                                sbk.write(line)          
                except:
                    pass
                
        