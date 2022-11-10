from distutils.log import info


pdbs = []
with open("/home/chengyj/kinase_work/dataset/Bridged_ring/pdb_rings/bridged.csv") as infos:
    for line in infos:
        info = line[0:-1].split("[")[1].replace("'","").split(",, ")
        for item in info:
            if item[0:4] not in pdbs:
                pdbs.append(item[0:4])
            else:
                pass
print(len(pdbs))

target = []
with open("/home/chengyj/kinase_work/dataset/Bridged_ring/pdb_rings/target.infos") as infos:
    for line in infos:
        if line.split("|")[1] in target:
            pass
        else:
            target.append(line.split("|")[1])

print(len(target))