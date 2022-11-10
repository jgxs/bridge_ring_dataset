from rdkit import Chem
from rdkit.Chem import Descriptors
patt = Chem.MolFromSmarts("[x3]")
def check_bridgedrings(mol_check,):
    hits = [item[0] for item in mol_check.GetSubstructMatches(patt)]
    hits_number = len(hits)
    for atom in hits:
        for id_atm in mol_check.GetAtomWithIdx(atom).GetNeighbors():
            if id_atm.GetIdx() in hits:
                hits_number -= 1
                break
    return hits_number > 0


patt_spiro = Chem.MolFromSmarts("[x4]")
patt_p_bridge = Chem.MolFromSmarts("[x3&!a]")
with open("/home/chengyj/kinase_work/ligands_from_chembl/bridged_rings_smis.csv") as data:
    for line in data:
        smis = line.split()[1]
        mol = Chem.MolFromSmiles(smis)
        ri = mol.GetRingInfo()
        largest_ring_size = max((len(r) for r in ri.AtomRings()), default=0)
        if largest_ring_size <= 12 and Descriptors.MolWt(mol) <= 1200:
            if mol.HasSubstructMatch(patt_spiro):
                with open("spiro.csv","a") as spiro:
                    spiro.write(line)
            elif check_bridgedrings(mol):
                with open("bridged.csv", "a") as bridge:
                    bridge.write(line)
            else:
                pass
        else:
            pass


            
