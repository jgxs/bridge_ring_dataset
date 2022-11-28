from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, rdFMCS
from rdkit.Chem import rdRGroupDecomposition as rdRGD
from tqdm import tqdm
from rdkit.Geometry.rdGeometry import Point3D
from pathlib import Path
# take off warnings
# 


def rename_atom(atom_to_rename, pdbinfo_ref):
    # sometime,
    name = str(atom_to_rename.GetSymbol()) + str(atom_to_rename.GetIdx() + 1)
    if len(name) == 2:
        name = " " + name + " "
    elif len(name) == 3:
        name = " " + name
    elif len(name) == 4:
        pass
    pdbinfo = atom_to_rename.GetPDBResidueInfo()
    if pdbinfo:
        Chem.AtomMonomerInfo.SetName(pdbinfo, name)
        Chem.AtomPDBResidueInfo.SetAltLoc(pdbinfo, " ")
        Chem.AtomPDBResidueInfo.SetResidueNumber(pdbinfo, 1)
    else:
        Chem.Atom.SetMonomerInfo(atom_to_rename, pdbinfo_ref)
        pdbinfo = atom_to_rename.GetPDBResidueInfo()
        Chem.AtomPDBResidueInfo.SetAltLoc(pdbinfo, " ")
        Chem.AtomMonomerInfo.SetName(pdbinfo, name)
    return name


def check_altLoc(altLoc_symbol, lines):
    # sometimes the altLoc in ligands line are not "A"," " or "1"
    # which leads to the pdbparser will ignore them
    # this function will select the ligand lines with the same altLoc.
    if altLoc_symbol == "A" or altLoc_symbol == " " or altLoc_symbol == "1":
        return "".join(lines)
    else:
        lig_Block = ""
        for line in lines:
            if line[16] == altLoc_symbol:
                lig_Block += line
        return lig_Block


def check_triazole_pdb(mol_pdb, smi):
    # sometimes the PDBparser will connect atoms uncorrectly
    # this function will remove the extra bond of the atom of "CN1N=NC=C1"
    mol_check = Chem.MolFromSmiles(smi)
    triazole = Chem.MolFromSmiles("CN1N=NC=C1")
    match_triazole = mol_check.GetSubstructMatch(triazole)
    if not match_triazole:
        return mol_pdb
    else:
        conf = mol_pdb.GetConformer()
        mol_edited = Chem.EditableMol(mol_pdb)
        dist_Idx = []
        for atom in mol_pdb.GetAtoms():
            atom_symbol = atom.GetSymbol()
            neighors = atom.GetNeighbors()
            if atom_symbol == "N" and len(neighors) > 3:
                N_idx = atom.GetIdx()
                dist_Idx = [
                    (
                        Point3D.Distance(
                            conf.GetAtomPosition(item.GetIdx()),
                            conf.GetAtomPosition(N_idx),
                        ),
                        item.GetIdx(),
                    )
                    for item in neighors
                ]
        if len(dist_Idx) > 0:
            dist_Idx.sort()
            mol_edited.RemoveBond(N_idx, dist_Idx[-1][1])
            return mol_edited.GetMol()
        else:
            return mol_pdb


def getconf_from_missingatom_pdb(refmol, inpmol):

    AllChem.EmbedMolecule(inpmol, maxAttempts=50)
    inpmol = Chem.RemoveAllHs(inpmol)

    mcs = rdFMCS.FindMCS(
        [inpmol, refmol],
        timeout=3,
        completeRingsOnly=True,
        bondCompare=rdFMCS.BondCompare.CompareAny,
    )
    if not mcs.queryMol:
        raise ValueError("there are too many missing atoms")
    bonded_conf = refmol.GetConformer()
    conf_res = inpmol.GetConformer()
    inpmol_prop = Chem.rdForceFieldHelpers.MMFFGetMoleculeProperties(inpmol)
    if inpmol_prop:
        ff_mcs = Chem.rdForceFieldHelpers.MMFFGetMoleculeForceField(inpmol, inpmol_prop)
    else:
        ff_mcs = Chem.rdForceFieldHelpers.UFFGetMoleculeForceField(inpmol)

    for i, j in zip(
        refmol.GetSubstructMatch(mcs.queryMol), inpmol.GetSubstructMatch(mcs.queryMol)
    ):
        ff_mcs.AddFixedPoint(j)
        conf_res.SetAtomPosition(j, bonded_conf.GetAtomPosition(i))

    for i in range(10):
        try:
            ff_mcs.Minimize()
            pass
        except:
            pass
    Chem.MolToPDBBlock(inpmol)
    return inpmol


def exctract_ligand_from_pdb(pdb_path, lig_id, smiles, outfile):
    # 抽取pdbfile中的ligand部分,获得rdkit可读字符串形式的pdbBlock
    # 并加氢重命名后以pdb格式储存在outfile中
    mol_from_smi = Chem.MolFromSmiles(smiles)
    num = mol_from_smi.GetNumAtoms()
    # print(num)
    ligand_lines = []
    with open(pdb_path) as ent:
        for line in ent:
            if line[0:6] == "HETATM" and line[17:20] == lig_id:
                ligand_lines.append(line[0:78] + "  \n")
    # print("".join(ligand_lines))
    ligand_start = 0
    ligand_end = len(ligand_lines)
    for i in range(1, len(ligand_lines)):
        if ligand_lines[i][17:26] != ligand_lines[i - 1][17:26]:
            # print(ligand_lines[i][17:26])
            ligand_end = i
            if ligand_end - ligand_start >= num:
                break
            else:
                ligand_start = i
                ligand_end = len(ligand_lines)
        else:
            pass
    # print(ligand_end)
    # print(ligand_start)
    altLoc = ligand_lines[ligand_start][16]
    if ligand_end - ligand_start <= 2 * num:
        lig_Block = check_altLoc(altLoc, ligand_lines[ligand_start:ligand_end])
    else:
        lig_Block = check_altLoc(
            altLoc, ligand_lines[ligand_start : ligand_start + num]
        )

    # print(lig_Block)

    # 加氢和重命名原子
    if altLoc == "A" or altLoc == " " or altLoc == "1":
        lig_pdb = Chem.MolFromPDBBlock(lig_Block)
    else:
        lig_pdb = Chem.MolFromPDBBlock(lig_Block, flavor=1)
    # print(lig_pdb.GetNumAtoms())
    if lig_pdb.GetNumAtoms() == 0:
        lig_pdb = Chem.MolFromPDBBlock(lig_Block, flavor=1)
        # print(lig_pdb.GetNumAtoms())

    lig_pdb = check_triazole_pdb(lig_pdb, smiles)

    if "H" in smiles:
        mol_from_smi = Chem.RemoveAllHs(mol_from_smi)
    # print("test")
    if lig_pdb.GetNumAtoms() >= mol_from_smi.GetNumAtoms():
        lig_pdb_refined = AllChem.AssignBondOrdersFromTemplate(mol_from_smi, lig_pdb)
    else:
        lig_pdb_refined = getconf_from_missingatom_pdb(lig_pdb, mol_from_smi)
    lig_pdb_refined.GetConformer()

    lig_pdb_H = Chem.AddHs(lig_pdb_refined, addCoords=True)

    atom_ref = lig_pdb.GetAtomWithIdx(1).GetPDBResidueInfo()
    Chem.AtomPDBResidueInfo.SetAltLoc(atom_ref, " ")
    Chem.AtomPDBResidueInfo.SetResidueNumber(atom_ref, 1)
    for atom in lig_pdb_H.GetAtoms():
        rename_atom(atom, atom_ref)
    Chem.MolToPDBFile(lig_pdb_H, outfile)

    return lig_pdb_H

def get_connecter(mol,patt):
    # in fact, this function is able to replace two function above.
    # however, i realise it too later... 
    match_patt = mol.GetSubstructMatches(patt)
    for match_set in match_patt:
        core_connecter = []
        side_connecter = []
        for idx in match_set:
            atom = mol.GetAtomWithIdx(idx)
            neighbors = {neighbor.GetIdx() for neighbor in atom.GetNeighbors()}
            if neighbors & set(match_set) == neighbors:
                if len(list(neighbors)) == 1:
                    side_connecter.append(idx)
                    core_connecter.append(list(neighbors)[0])
                else:
                    pass
            else:
                side_connecter.append(idx)
                for i in list(neighbors):
                    if i in match_set:
                        core_connecter.append(i)
        if len(side_connecter) == 2:
            return [core_connecter, side_connecter,list(set(match_set)-set(side_connecter))]
        else:
            pass
    return False


def phe2bch_with_ligfpbd(lig,output):
    lig_noH = Chem.RemoveHs(lig)
    conf = lig_noH.GetConformer()
    lig_noH_atoms = lig_noH.GetNumAtoms()

    patt = Chem.MolFromSmarts("*!:c1cc(!:*)ccc1")
    connect_atoms = get_connecter(lig_noH, patt)

    re_pat = Chem.MolFromSmiles("C12CCCC(C2)C1")
    AllChem.EmbedMolecule(re_pat)
    mid_mol = Chem.CombineMols(lig_noH, re_pat)

    mid_edited = Chem.EditableMol(mid_mol)
    mid_edited.AddBond(
        connect_atoms[1][0], lig_noH_atoms, order=Chem.rdchem.BondType.SINGLE
    )
    mid_edited.AddBond(
        connect_atoms[1][1], lig_noH_atoms + 4, order=Chem.rdchem.BondType.SINGLE
    )
    try:
        mid_edited.RemoveBond(connect_atoms[1][0],connect_atoms[0][0])
        mid_edited.RemoveBond(connect_atoms[1][1],connect_atoms[0][1])
    except:
        mid_edited.RemoveBond(connect_atoms[1][0],connect_atoms[0][1])
        mid_edited.RemoveBond(connect_atoms[1][1],connect_atoms[0][0])
    phenyl_atoms = connect_atoms[2]
    phenyl_atoms.sort(reverse=True)
    for i in phenyl_atoms:
        mid_edited.RemoveAtom(i)
    BCHep_mol = mid_edited.GetMol()
    Chem.SanitizeMol(BCHep_mol)
    conf_gen = BCHep_mol.GetConformer()
    mol_prop = Chem.rdForceFieldHelpers.MMFFGetMoleculeProperties(BCHep_mol)
    if mol_prop:
        ff_mcs = Chem.rdForceFieldHelpers.MMFFGetMoleculeForceField(BCHep_mol, mol_prop)
    else:
        ff_mcs = Chem.rdForceFieldHelpers.UFFGetMoleculeForceField(BCHep_mol)
    for atom in BCHep_mol.GetAtoms():
        atom_idx = atom.GetIdx()
        if atom_idx < lig_noH_atoms - 6:
            ff_mcs.AddFixedPoint(atom_idx)
        elif atom_idx == lig_noH_atoms-6:
            ff_mcs.AddFixedPoint(atom_idx)
            conf_gen.SetAtomPosition(atom_idx, conf.GetAtomPosition(list(connect_atoms[0])[0]))
        elif atom_idx == lig_noH_atoms-2:
            ff_mcs.AddFixedPoint(atom_idx)
            conf_gen.SetAtomPosition(atom_idx, conf.GetAtomPosition(list(connect_atoms[0])[1]))
        else:
            pass
    for i in range(100):
        try:
            ff_mcs.Minimize()
        except:
            pass
    BCHep_mol_H = Chem.AddHs(BCHep_mol, addCoords=True)
    Chem.MolToPDBFile(BCHep_mol_H, output)
    return BCHep_mol_H


def sin_work(inp,ligands_smi):
    if Path("sinwork.err").exists():
        pass
    else:
        with open("sinwork.err","w") as log:
            log.write("test begin\n")
    key = inp.split("_")[0]
    pdb_id = inp.split("_")[1]
    lig_smi = ligands_smi[key][0]
    pdb_path = f"/home/chengyj/kinase_work/dataset/Bridged_ring/PDB_rings/PHE2BCH_pairs/pdb_dataset/pdb/pdb{pdb_id}.ent"
    try:
        lig_Block = exctract_ligand_from_pdb(pdb_path, key, lig_smi, f"{key}_{pdb_id}_phe.pdb")
    except Exception as ex:
        with open("sinwork.err", "a") as log:
            log.write(f"PDBload error of {key}_{pdb_id}: {ex}" + "\n")
        return False
    try:
        BCH_lig_H = phe2bch_with_ligfpbd(lig_Block, f"{key}_{pdb_id}_bch.pdb")
    except Exception as ex:
        with open("sinwork.err", "a") as log:
            log.write(f"construct error of {key}_{pdb_id}: {ex}" + "\n")
        return False
    return BCH_lig_H


if __name__ == "__main__":
    RDLogger.DisableLog("rdApp.*")
    ligands_smi = {}
    with open(
        "/home/chengyj/kinase_work/dataset/Bridged_ring/PDB_rings/PHE2BCH_pairs/scripts/AaaaA_only5.csv"
    ) as AaaaA:
        for line in AaaaA:
            info = line.split()
            ligands_smi[info[1]] = [info[0], False]
    # print(1)
    with open(
        "/home/chengyj/kinase_work/dataset/Bridged_ring/PDB_rings/PHE2BCH_pairs/pdb_dataset/cc-to-pdb.tdd"
    ) as ligs:
        for line in ligs:
            info = line.split()
            if info[0] in ligands_smi.keys():
                ligands_smi[info[0]].append(info[1:])
                ligands_smi[info[0]][1] = True
    # print(2)

    print(sin_work("11B_2pj5", ligands_smi))
    
    for key in tqdm(ligands_smi):
        if ligands_smi[key][1]:
            for pdbid in ligands_smi[key][2]:
                # print(f"{key}_{pdbid}")
                sin_work(f"{key}_{pdbid}", ligands_smi)
                """
                with open("phe2bch.log", "a") as log:
                    log.write(f"{pdbid}\n")
                pdb_file_path = pdb_dir + f"pdb{pdbid}.ent"
                sin_work(f"{key}_{pdbid}", ligands_smi)
                try:
                    lig_Block = exctract_ligand_from_pdb(
                        pdb_file_path,
                        key,
                        ligands_smi[key][0],
                        f"{key}_{pdbid}_phe.pdb",
                    )
                except Exception as ex:
                    with open("phe2bch.err", "a") as log:
                        log.write(f"PDBload error of {key}_{pdbid}: {ex}" + "\n")

                #try:
                print(f"{key}_{pdbid}")
                test = phe2bch_topdb(
                        ligands_smi[key][0], lig_Block, f"{key}_{pdbid}_bch.pdb"
                    )
                print("finish")
                #except Exception as ex:
                #    with open("phe2bch.err", "a") as log:
                #        log.write(f"Error of {key}_{pdbid}: {ex}" + "\n")
                """ 