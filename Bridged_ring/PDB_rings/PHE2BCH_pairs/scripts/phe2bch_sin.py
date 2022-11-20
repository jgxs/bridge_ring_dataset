from Phe2BCH import *


ligands_smi = {}
with open(
    "/home/chengyj/kinase_work/dataset/Bridged_ring/PDB_rings/PHE2BCH_pairs/lig_menu/AaaaA_only3.csv"
) as AaaaA:
    for line in AaaaA:
        info = line.split()
        ligands_smi[info[1]] = [info[0]]
# inp = "0FS_3vc4" segmentation fault
inp = "RFQ_4d3u"
key = inp.split("_")[0]
pdb_id = inp.split("_")[1]
lig_smi = ligands_smi[key][0]
# print(lig_smi)
pdb_path = f"/home/chengyj/kinase_work/dataset/Bridged_ring/PDB_rings/PHE2BCH_pairs/pdb_dataset/pdb/pdb{pdb_id}.ent"


lig_Block = exctract_ligand_from_pdb(pdb_path, key, lig_smi, f"{key}_{pdb_id}_sin.pdb")
# print("test")
test = phe2bch_topdb(lig_smi, lig_Block, f"{key}_{pdb_id}_bch.pdb")
