from Phe2BCH import *


ligands_smi = {}
with open("/home/chengyj/kinase_work/dataset/Bridged_ring/PDB_rings/PHE2BCH_pairs/lig_menu/AaaaA_only3.csv") as AaaaA:
    for line in AaaaA:
        info = line.split()
        ligands_smi[info[1]] = [info[0]]
key='PIX'
pdb_id='3lj6'
lig_smi=ligands_smi[key][0]
# print(lig_smi)
pdb_path = f"/home/chengyj/kinase_work/dataset/Bridged_ring/PDB_rings/PHE2BCH_pairs/pdb_dataset/pdb/pdb{pdb_id}.ent"
lig_Block = exctract_ligand_from_pdb(pdb_path,key,lig_smi,"test_sin.pdb")
# print(lig_Block)
# print(lig_Block)
with open(f"{key}_{pdb_id}_phe.pdb","w") as phe:
    phe.write(lig_Block)
test = phe2bch_topdb(lig_smi,lig_Block,f"{key}_{pdb_id}_bch.pdb")



