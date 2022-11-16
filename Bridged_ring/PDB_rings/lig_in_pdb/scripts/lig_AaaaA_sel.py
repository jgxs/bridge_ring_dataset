from tqdm import tqdm
from rdkit import Chem, RDLogger
# take off warnings
RDLogger.DisableLog('rdApp.*')

patt_sel = Chem.MolFromSmarts("*!:c1cc(!:*)ccc1")
patt_not1 = Chem.MolFromSmarts("*c1c(*)cccc1")
patt_not2 = Chem.MolFromSmarts("*c1ccc(*)cc1")
patt_not3 = Chem.MolFromSmarts("*ccc(*)cc(*)c")
patt_not4 = Chem.MolFromSmarts("*ccc(*)c(*)cc")

ligand_r6 = []
with open("/home/chengyj/kinase_work/dataset/Bridged_ring/PDB_rings/lig_in_pdb/pdb_dataset/Components-smiles-stereo-oe.smi") as inputfile:
    for line in tqdm(inputfile):
        try:
            mol = Chem.MolFromSmiles(line.split()[0])
            if mol.GetSubstructMatches(patt_sel) and not mol.GetSubstructMatches(patt_not1) and not mol.GetSubstructMatches(patt_not2) and not mol.GetSubstructMatches(patt_not3)and not mol.GetSubstructMatches(patt_not4):
                ligand_r6.append(line)
        except:
            pass

with open("lig_menu/AaaaA_only4.csv","w") as output:
    for line in ligand_r6:
        output.write(line)