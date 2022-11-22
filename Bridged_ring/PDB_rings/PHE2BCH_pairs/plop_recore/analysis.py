from pathlib import Path

rescore_dir = Path("/pubhome/yjcheng02/bridge_ring_dataset/Bridged_ring/PDB_rings/PHE2BCH_pairs/lig_bch_refine")
all_res = list(rescore_dir.glob("**/res.out"))

with open("rescore.log","w") as log:
    log.write("analysis starting \n")
with open("rescore.err","w") as log:
    log.write("analysis starting \n")
with open("rescore.better","w") as log:
    log.write("")

for res in all_res:
    lig = res.parts[-3]
    pdb = res.parts[-2]
    with res.open() as infom:
        info = infom.read().split("\n")
    phe_score = float(info[0].split()[1])
    bch_score = float(info[1].split()[1])
    if phe_score >= 0.0 or bch_score >= 0.0:
        with open("rescore.err","a") as err:
            err.write(f"Error of {lig}_{pdb}\n")
    else:
        with open("rescore.log","a") as log:
            log.write(f"{lig} {pdb} {phe_score: >12.4f} {bch_score: >12.4f} {bch_score-phe_score: >12.4f}\n")
        if bch_score-phe_score < 0:
            with open("rescore.better","a") as log:
                log.write(f"{lig} {pdb} {phe_score: >12.4f} {bch_score: >12.4f} {bch_score-phe_score: >12.4f}\n")
     

