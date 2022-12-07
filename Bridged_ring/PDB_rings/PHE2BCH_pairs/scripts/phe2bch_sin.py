from Phe2BCH import *

RDLogger.DisableLog("rdApp.*")
ligands_smi = {}
with open(
    "/home/chengyj/kinase_work/dataset/Bridged_ring/PDB_rings/PHE2BCH_pairs/lig_menu/AaaaA_small_rings.csv"
) as AaaaA:
    for line in AaaaA:
        info = line.split()
        ligands_smi[info[1]] = [info[0]]

# inp = "0FS_3vc4" segmentation fault
# print("test")
inps = ["0OO_4eh7","0RZ_4f1q","07Q_3tv7","B74_6xic","0VO_4g1y","9F5_5o0d","0BU_3upf","EDR_2ack","EDR_1ax9","IMA_1lqe","IMA_1lpg","NAF_2h9y","NAF_1amn","MTY_1biq"]
inps = ["MTY_1biq"]
inps = []
with open("menu") as info:
    for line in info:
        inps.append(line[0:-1])
for inp in inps:
    sin_work(inp, ligands_smi)
# sin_work("11B_2pj5", ligands_smi)