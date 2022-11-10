import psycopg2
import psycopg2.extras
# https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/schema_documentation.html
# key value table

from rdkit import Chem
from rdkit.Chem import Descriptors

patt = Chem.MolFromSmarts("[x3]")
patt_spiro = Chem.MolFromSmarts("[x4]")
patt_p_bridge = Chem.MolFromSmarts("[x3&!a]")

couts = [0]*5
sum_all = 0
bigrings = 1
sprio = 2
bridge =3
spk = 4
print(couts)
with open("/home/chengyj/kinase_work/Bridged_ring/chembl_rings/molregno_smis.csv") as mols:
    for mol in mols:
        couts[sum_all] += 1
        smis = mol.split()[1]
        try:
            mol = Chem.MolFromSmiles(smis)
            ri = mol.GetRingInfo()
            largest_ring_size = max((len(r) for r in ri.AtomRings()), default=0)
            if largest_ring_size > 12 and Descriptors.MolWt(mol) <= 1200:
                couts[bigrings] += 1
            else: 
                if mol.HasSubstructMatch(patt_spiro):
                    couts[sprio] += 1
                if mol.HasSubstructMatch(patt_p_bridge):
                    couts[bridge] += 1        
                if mol.HasSubstructMatch(patt_spiro) and mol.HasSubstructMatch(patt_p_bridge):
                    couts[spk] += 1
        except:
            pass
print(couts)
        