from pathlib import Path

def dG_calculate(line):
    if "0.00" not in line and "nan" not in line and "NaN" not in line:
        data = line.split()
        cmx = float(data[1])
        rec = float(data[6])
        lig = float(data[10])
        return cmx - rec -lig
    else:
        return 0.00

rescore_dir = Path(".")
rescore_name = rescore_dir.absolute().name
print(rescore_name)
all_res = list(rescore_dir.glob("**/res.out"))
res_dir = "/pubhome/yjcheng02/bridge_ring_dataset/Bridged_ring/PDB_rings/PHE2BCH_pairs/charged_res/"

Herr = []
Eerr = []
better = []
normol = []
for res in all_res:
    lig = res.parts[-3]
    pdb = res.parts[-2]
    with res.open() as infom:
        info = infom.read().split("\n")
    phe_score = dG_calculate(info[0])
    bch_score = dG_calculate(info[1])
    if abs(phe_score) > 1000:
        Herr.append(f"{lig} {pdb} {phe_score: >24.4f} {bch_score: >24.4f} {bch_score-phe_score: >12.4f}\n")
    elif phe_score == 0.0 or bch_score == 0.0:
        Eerr.append(f"Error of molecule {lig}_{pdb} \n")
    else:
        normol.append(f"{lig}_{pdb} {phe_score: >12.4f} {bch_score: >12.4f} {bch_score-phe_score: >12.4f}\n")
        if phe_score < 0.0 and bch_score < phe_score:
            better.append(f"{lig}_{pdb} {phe_score: >12.4f} {bch_score: >12.4f} {bch_score-phe_score: >12.4f}\n")

outlog = {"log":normol,"Herr":Herr,"Eerr":Eerr,"better":better}
print(f"result is in {res_dir}")
for key in outlog:
    with open(f"{res_dir}{rescore_name}.{key}","w") as log:
        for line in outlog[key]:
            log.write(line)

