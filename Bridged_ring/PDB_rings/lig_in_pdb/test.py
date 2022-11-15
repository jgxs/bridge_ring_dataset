from Phe2BCH import *


ligands_smi = {}
with open("/home/chengyj/kinase_work/dataset/Bridged_ring/PDB_rings/lig_in_pdb/lig_menu/AaaaA_only3.csv") as AaaaA:
    for line in AaaaA:
        info = line.split()
        ligands_smi[info[1]] = [info[0]]
key='5N0'
pdb_id='7ril'
lig_smi=ligands_smi[key][0]
# print(lig_smi)
pdb_path = f"/home/chengyj/kinase_work/dataset/Bridged_ring/PDB_rings/lig_in_pdb/pdb_dataset/pdb/pdb{pdb_id}.ent"
lig_Block = pdb_2_lig_block(pdb_path,key,lig_smi)
# print(lig_Block)
# print(lig_Block)
test = phe2bch_topdb(lig_smi,lig_Block,f"{key}_{pdb_id}.pdb")



