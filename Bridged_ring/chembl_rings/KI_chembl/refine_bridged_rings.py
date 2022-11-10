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
with open("ki_chembl.smi") as data:
    for line in data:
        smis = line.split()[1]
        mol = Chem.MolFromSmiles(smis)
        ri = mol.GetRingInfo()
        largest_ring_size = max((len(r) for r in ri.AtomRings()), default=0)
        if largest_ring_size > 12 and Descriptors.MolWt(mol) <= 1200:
            with open("bigrings_ki.csv","a") as bigring:
                bigring.write(line)
        else: 
            if mol.HasSubstructMatch(patt_spiro):
                with open("spiro_ki.csv","a") as spiro:
                    spiro.write(line)
            if check_bridgedrings(mol):
                with open("bridged_ki.csv", "a") as bridge:
                    bridge.write(line)
            else:
                pass
            
