from pathlib import Path
from tqdm import tqdm
from pdb_parser import Ligand


def ligand2chain_id(path):
    with open(path) as info:
        for line in info:
            return line[21]


def split_target_chain(pdbfile, chain_id, target_path):
    pdb_block = []
    with open(pdbfile) as pdblines:
        for line in pdblines:
            if line[0:4] == "ATOM" and line[21] == chain_id:
                pdb_block.append(line)
    with open(target_path, "w") as output:
        for item in pdb_block:
            output.write(item)
        output.write("TER\n")

with open("check.log", "w") as log:
    log.write("check staring\n")

all_ligs_str = "/home/chengyj/kinase_work/dataset/Bridged_ring/PDB_rings/PHE2BCH_pairs/ligs_charged/"
all_ligs_dir = Path(all_ligs_str)

all_ligs_pdb = list(all_ligs_dir.glob("**/*bch.pdb"))
lig_dict = {}
for item in all_ligs_pdb:
    ligid, pdbid, type = item.name.replace(".pdb", "").split("_")
    phe_path = all_ligs_dir/f"{ligid}_{pdbid}_phe.pdb"
    if phe_path.exists():
        if ligid not in lig_dict:
            lig_dict[ligid] = []
            lig_dict[ligid].append(pdbid)
        else:
            lig_dict[ligid].append(pdbid)
    else:
        with open("check.log", "w") as log:
            log.write(f"{ligid}_{pdbid}_phe.pdb is missing\n")
        pass

lig_complex_dir = str(
    "/home/chengyj/kinase_work/dataset/Bridged_ring/PDB_rings/PHE2BCH_pairs/lig_bch_charged"
)
try:
    Path(lig_complex_dir).mkdir()
except:
    pass

for key in tqdm(lig_dict):
    for value in lig_dict[key]:
        # check the ligand is right
        ligand_stat = True
        lig_phe_bch = {}
        for group in ["bch", "phe"]:
            lig_complex = Path(
                f"/home/chengyj/kinase_work/dataset/Bridged_ring/PDB_rings/PHE2BCH_pairs/ligs_charged/{key}_{value}_{group}.pdb"
            )
            lig_phe_bch[group] = Ligand(lig_complex)

        for atom_bch_idx in lig_phe_bch["bch"].atoms:
            min_dist = 100
            atom_bch = lig_phe_bch["bch"].atoms[atom_bch_idx]
            for atom_phe_idx in lig_phe_bch["phe"].atoms:
                atom_phe = lig_phe_bch["phe"].atoms[atom_phe_idx]
                local_dist = atom_bch.get_dist(atom_phe)
                min_dist = local_dist if local_dist < min_dist else min_dist
            if min_dist > 20:
                ligand_stat = False
                with open("check.log", "a") as output:
                    output.write(f"{key}_{value}_bch.pdb is wrong\n")
                break

        if ligand_stat:
            ligs = Path(f"{lig_complex_dir}/{key}")
            try:
                ligs.mkdir()
            except:
                pass
            complex_dir = Path(f"{lig_complex_dir}/{key}/{value}")
            try:
                complex_dir.mkdir()
            except:
                pass
            lig_phe_bch = {}
            for group in ["bch", "phe"]:
                lig_complex = Path(
                    f"{lig_complex_dir}/{key}/{value}/{key}_{value}_{group}.pdb"
                )
                lig_complex.symlink_to(
                    f"/home/chengyj/kinase_work/dataset/Bridged_ring/PDB_rings/PHE2BCH_pairs/ligs_charged/{key}_{value}_{group}.pdb"
                )

            pdbent = f"/home/chengyj/kinase_work/dataset/Bridged_ring/PDB_rings/PHE2BCH_pairs/pdb_dataset/pdb/pdb{value}.ent"
            # pdbpath = Path(f"{lig_complex_dir}/{key}/{value}/{value}.pdb")
            # pdbpath.symlink_to(pdbent)

            lig_phe = f"{lig_complex_dir}/{key}/{value}/{key}_{value}_phe.pdb"
            targer_chain_id = ligand2chain_id(lig_phe)
            split_target_chain(
                pdbent,
                targer_chain_id,
                f"{lig_complex_dir}/{key}/{value}/{value}_sin.pdb",
            )
