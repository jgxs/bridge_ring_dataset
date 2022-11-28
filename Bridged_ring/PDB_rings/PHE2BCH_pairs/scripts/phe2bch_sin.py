from Phe2BCH import *


ligands_smi = {}
with open(
    "/home/chengyj/kinase_work/dataset/Bridged_ring/PDB_rings/PHE2BCH_pairs/lig_menu/AaaaA_only3.csv"
) as AaaaA:
    for line in AaaaA:
        info = line.split()
        ligands_smi[info[1]] = [info[0]]

# inp = "0FS_3vc4" segmentation fault
# print("test")
inp = "0OO_4eh7"
# inp = "0RZ_4f1q"
inp = "07Q_3tv7"
# inp = "B74_6xic"

sin_work(inp, ligands_smi)
sin_work("11B_2pj5", ligands_smi)