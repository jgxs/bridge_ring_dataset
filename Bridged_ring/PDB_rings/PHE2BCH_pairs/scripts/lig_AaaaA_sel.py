from tqdm import tqdm
from rdkit import Chem, RDLogger

# take off warnings
RDLogger.DisableLog("rdApp.*")

patt_sel = Chem.MolFromSmarts("*!:c1cc(!:*)ccc1")

def get_connecter(mol,patt):
    # in fact, this function is able to replace two function above.
    # however, i realise it too later... 
    match_patt = mol.GetSubstructMatches(patt)
    for match_set in match_patt:
        core_connecter = []
        side_connecter = []
        for idx in match_set:
            atom = mol.GetAtomWithIdx(idx)
            neighbors = {neighbor.GetIdx() for neighbor in atom.GetNeighbors()}
            if neighbors & set(match_set) == neighbors:
                if len(list(neighbors)) == 1:
                    side_connecter.append(idx)
                    core_connecter.append(list(neighbors)[0])
                else:
                    pass
            else:
                side_connecter.append(idx)
                for i in list(neighbors):
                    if i in match_set:
                        core_connecter.append(i)
        # print(side_connecter)
        if len(side_connecter) == 2:
            return [core_connecter, side_connecter,list(set(match_set)-set(side_connecter))]
        else:
            pass
    return False

ligand_r6 = []

with open(
    "/home/chengyj/kinase_work/dataset/Bridged_ring/PDB_rings/PHE2BCH_pairs/lig_menu/Components-smiles-stereo-oe.smi"
) as inputfile:
    for line in tqdm(inputfile):
        try:
            mol = Chem.MolFromSmiles(line.split()[0])
            ri = mol.GetRingInfo()
            largest_ring_size = max((len(r) for r in ri.AtomRings()), default=0)
            if (
                mol.GetSubstructMatches(patt_sel) 
                and get_connecter(mol,patt_sel)
                and largest_ring_size < 14
            ):
                ligand_r6.append(line)
        except:
            pass

with open("AaaaA_only5.csv", "w") as output:
    for line in ligand_r6:
        output.write(line)
