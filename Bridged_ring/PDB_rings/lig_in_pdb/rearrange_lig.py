from pathlib import Path

all_ligs_dir = Path("/home/chengyj/kinase_work/dataset/Bridged_ring/PDB_rings/lig_in_pdb/bak/AaaaA4_ref/")

all_ligs_pdb = list(all_ligs_dir.glob('**/*bch.pdb'))
lig_dict = {}
for item in all_ligs_pdb:
    ligid,pdbid,type = item.parts[-1].replace(".pdb","").split("_")
    if ligid not in lig_dict:
        lig_dict[ligid] = []
        lig_dict[ligid].append(pdbid)
    else:
        lig_dict[ligid].append(pdbid)
lig_complex_dir=str("/home/chengyj/kinase_work/dataset/Bridged_ring/PDB_rings/lig_in_pdb/lig_bch_complex")
for key in lig_dict:
    ligs = Path(f"{lig_complex_dir}/{key}")
    ligs.mkdir()
    for value in lig_dict[key]:
        complex = Path(f"{lig_complex_dir}/{key}/{value}")
        complex.mkdir()
        for group in ["bch","phe"]:
            lig_complex = Path(f"{lig_complex_dir}/{key}/{value}/{key}_{value}_{group}.pdb")
            lig_complex.symlink_to(f"/home/chengyj/kinase_work/dataset/Bridged_ring/PDB_rings/lig_in_pdb/bak/AaaaA4_ref/{key}_{value}_{group}.pdb")
        pdbent = f"/home/chengyj/kinase_work/dataset/Bridged_ring/PDB_rings/lig_in_pdb/pdb_dataset/pdb/pdb{value}.ent"
        pdbpath = Path(f"{lig_complex_dir}/{key}/{value}/{value}.pdb")
        pdbpath.symlink_to(pdbent)

    
