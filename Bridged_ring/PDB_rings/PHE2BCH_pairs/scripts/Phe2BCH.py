from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, rdFMCS
from rdkit.Chem import rdRGroupDecomposition as rdRGD
from tqdm import tqdm
from rdkit.Geometry.rdGeometry import Point3D

# take off warnings
RDLogger.DisableLog("rdApp.*")


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


def get_connect_atom_in_core(orimol, sidechain_atom_idxes):
    # mol is the whole molecule, the sidechain_atom_idx is the a tuple of the atom idx of a sidechain;
    # this function will get the atom connected with the sidechain in the core of the orimol.
    connect_atoms_in_core = []
    for idx in sidechain_atom_idxes:
        for atom in orimol.GetAtomWithIdx(idx).GetNeighbors():
            if atom.GetIdx() not in sidechain_atom_idxes:
                connect_atoms_in_core.append(atom.GetIdx())
    return connect_atoms_in_core


def get_connect_atom_in_sidechain(orimol, core):
    # the core is the a tuple of the atom idx of the core;
    # this function will get the atom connected with the core in the sidechain of the orimol.
    rings = orimol.GetRingInfo()
    for r in rings.AtomRings():
        count = 0
        for idx in core:
            if idx in set(r):
                count += 1
        if count == len(core) - 2:
            core_atom_idx = set(r)
            connect_atom_in_sidechain = set(core) - core_atom_idx
    return list(connect_atom_in_sidechain)

def get_connecter(mol,patt):
    # in fact, this function is able to replace two function above.
    # however, i realise it too later... 
    core_connecter = []
    side_connecter = []
    match_patt = mol.GetSubstructMatches(patt)
    match_set = set(match_patt[0])
    for idx in match_patt[0]:
        atom = mol.GetAtomWithIdx(idx)
        neighbors = {neighbor.GetIdx() for neighbor in atom.GetNeighbors()}
        if neighbors & match_set == neighbors:
            pass
        else:
            side_connecter.append(idx)
            for i in list(neighbors):
                if i in match_set:
                    core_connecter.append(i)
    return [set(core_connecter),set(side_connecter)]

def phe2bch_with_smiles(smi):
    # 用于将给定配体中的*间位苯环骨架*替换为*螺旋桨烷骨架*的脚本
    # 用smi指定ligand的结构

    # 读入分子
    mol = Chem.MolFromSmiles(smi)
    # 识别间位取代部分，* 匹配任意原子；'!:' 非芳香键； *!:ccc(!:*)ccc 含间位取代的苯环
    patt = Chem.MolFromSmarts("*!:c1cc(!:*)ccc1")
    match = mol.GetSubstructMatches(patt)

    # 识别两个间位取代基与苯环的连接原子,识别需要被替换苯环原子
    connect_atoms = get_connect_atom_in_sidechain(mol, match[0])
    core = set(match[0]) - set(connect_atoms)

    # 记录第一部分的原子数目
    mol1_Atom_num = mol.GetNumAtoms()
    # 引入螺旋桨烷，并将其与原分子置于一处，使得其编号从·mol1_Atom_num+0·开始
    re_pat = Chem.MolFromSmiles("C12CCCC(C2)C1")
    mid_mol = Chem.CombineMols(mol, re_pat)
    # 开始成键
    mid_edited = Chem.EditableMol(mid_mol)
    # 螺旋桨烷烃的桥头碳分别为mol1_Atom_num和mol1_Atom_num+4
    mid_edited.AddBond(
        connect_atoms[0], mol1_Atom_num, order=Chem.rdchem.BondType.SINGLE
    )
    mid_edited.AddBond(
        connect_atoms[1], mol1_Atom_num + 4, order=Chem.rdchem.BondType.SINGLE
    )
    # 删除原来的苯环
    # 每删除一个原子后原子序号都会发生改变，所以需要倒序删除。
    phenyl_atoms = list(core)
    phenyl_atoms.sort(reverse=True)
    for i in phenyl_atoms:
        mid_edited.RemoveAtom(i)
    BCHep_mol = mid_edited.GetMol()
    # check mol
    try:
        Chem.SanitizeMol(BCHep_mol)
    except:
        raise ValueError("no match 6-5 rings")
    return BCHep_mol

def get_correct_match(mol_whole, mcs, arranged_atoms, core_mol, debug):
    match_to_use = False
    correct_matches = mol_whole.GetSubstructMatches(mcs.queryMol)
    if debug:
        print(f"begin {correct_matches}")
        print(arranged_atoms)
    if len(correct_matches) > 1:
        if debug:
            print(f"begin2 {correct_matches}")
        # arranged_atoms is not empty means it is the first bigger group to match
        # however, the len(correct_matches) > 1 means the bigger group is not 
        # big enough to match only one part of the refmol
        # In that situation, we need to select the correct part of the molecule.
        core_head, side_head = get_connecter(mol_whole,core_mol)
        for match_correct in correct_matches:
            #print(set(match_correct) & core_head)
            #print(set(match_correct) & side_head)
            if set(match_correct) & arranged_atoms == set() and set(match_correct) & core_head == set() and set(match_correct) & side_head != side_head:
                match_to_use = match_correct
    else:
        if debug:
            print(f"begin4 {correct_matches}")
        match_to_use = correct_matches[0]
    if debug:
        print(f"end {match_to_use}")
    if match_to_use:
        return match_to_use
    else:
        return correct_matches[0]


def getpdb(refmol, inpmol, pdbfile):

    bonded_conf = refmol.GetConformer()
    conf_res = inpmol.GetConformer()
    inpmol_prop = Chem.rdForceFieldHelpers.MMFFGetMoleculeProperties(inpmol)
    if inpmol_prop:
        ff_mcs = Chem.rdForceFieldHelpers.MMFFGetMoleculeForceField(inpmol, inpmol_prop)
    else:
        ff_mcs = Chem.rdForceFieldHelpers.UFFGetMoleculeForceField(inpmol)

    ref_arranged_atoms = set()
    inp_arranged_atoms = set()
    core = Chem.MolFromSmiles("C12CCCC(C2)C1")
    bch_core = Chem.MolFromSmarts("*C12CCCC(*)(C2)C1")
    phe_core = Chem.MolFromSmarts("*!:c1cc(!:*)ccc1")
    # maybe a better way while i do not want test it
    # bch_core_connect, bch_side_connect = get_connecter(inpmol,bch_core)
    # phe_core_connect, phe_side_connect = get_connecter(refmol,phe_core)

    side_chains, unmatched = rdRGD.RGroupDecompose([core], [inpmol])
    if side_chains[0]["R1"].GetNumAtoms()>side_chains[0]["R2"].GetNumAtoms():
        side_chains_reorder = [side_chains[0]["R1"],side_chains[0]["R2"]]  
    else:
        side_chains_reorder = [side_chains[0]["R2"],side_chains[0]["R1"]]
    
    for sidechain in side_chains_reorder:
        # Chem.SanitizeMol(sidechain)
        # print(ref_arranged_atoms)
        # print(inp_arranged_atoms)
        #for item in sidechain.GetAtoms():
        #    print(item.GetSymbol())
        sidechain_conf = sidechain.GetConformer()
        mcs1 = rdFMCS.FindMCS([refmol, sidechain], timeout=3, bondCompare=rdFMCS.BondCompare.CompareAny)
        mcs3 = rdFMCS.FindMCS([inpmol, sidechain], timeout=3, bondCompare=rdFMCS.BondCompare.CompareOrderExact)

        sidechain_ref = get_correct_match(refmol,mcs1,ref_arranged_atoms,phe_core,False)
        sidechain_gen = get_correct_match(inpmol,mcs3,inp_arranged_atoms,bch_core,False)
        
        atom_connect_sidechain_core_ref = get_connect_atom_in_core(refmol, sidechain_ref)[0]
        atom_connect_sidechain_core_gen = get_connect_atom_in_core(inpmol, sidechain_gen)[0]

        ff_mcs.AddFixedPoint(atom_connect_sidechain_core_gen)
        conf_res.SetAtomPosition(atom_connect_sidechain_core_gen,bonded_conf.GetAtomPosition(atom_connect_sidechain_core_ref),)
        
        mcs_ref2mid = rdFMCS.FindMCS([refmol, sidechain], timeout=3, bondCompare=rdFMCS.BondCompare.CompareAny)
        match_to_correct_ref = get_correct_match(refmol,mcs_ref2mid,ref_arranged_atoms,phe_core,False)
        # print(f"the result: {sidechain.GetSubstructMatches(mcs_ref2mid.queryMol)}")
        for i, j in zip(sidechain.GetSubstructMatch(mcs_ref2mid.queryMol),match_to_correct_ref,):
            sidechain_conf.SetAtomPosition(i, bonded_conf.GetAtomPosition(j))
            # print(f"set atom idx {i}")
            ref_arranged_atoms.add(j)
        # print(ref_arranged_atoms)
        mcs_mid2inp = rdFMCS.FindMCS(
            [inpmol, sidechain],
            timeout=3,
            bondCompare=rdFMCS.BondCompare.CompareOrderExact,
        )
        match_to_correct_inp = get_correct_match(inpmol,mcs_mid2inp,inp_arranged_atoms,bch_core,False)
        # print(match_to_correct_inp)
        # print(f"sidechain {sidechain.GetSubstructMatches(mcs_mid2inp.queryMol)}")
        for i, j in zip(
            sidechain.GetSubstructMatch(mcs_mid2inp.queryMol),
            match_to_correct_inp
        ):  
            # print(f"{i}      i")
            ff_mcs.AddFixedPoint(j)
            conf_res.SetAtomPosition(j,sidechain_conf.GetAtomPosition(i))
            inp_arranged_atoms.add(j)

    for i in range(50):
        try:
            ff_mcs.Minimize()
        except:
            pass

    inpmol = Chem.AddHs(inpmol, addCoords=True)
    Chem.MolToPDBFile(inpmol, pdbfile)
    return inpmol


def phe2bch_topdb(smi0, refpdb, name):
    mol = phe2bch_with_smiles(smi0)
    ## Totally can not understand
    # here must mol to smi to mol
    # otherwise, the Allchem.EmbedMolecule may be
    # Segmentation fault
    # the bug detail is the "0FS_3vc4"
    smi = Chem.MolToSmiles(mol)
    mol = Chem.MolFromSmiles(smi)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, maxAttempts=5000)
    # print("test")
    mol = Chem.RemoveAllHs(mol)

    structure_from_pdb = refpdb
    structure_from_pdb.GetConformer()
    try:
        structure_refine = Chem.RemoveAllHs(structure_from_pdb)
    except:
        structure_refine = Chem.RemoveHs(structure_from_pdb, implicitOnly=True)
    # some ligand can not remove all H
    structure_from_pdb.GetConformer()
    # print("test")
    return getpdb(structure_refine, mol, f"{name}")


if __name__ == "__main__":
    ligands_smi = {}
    with open(
        "/home/chengyj/kinase_work/dataset/Bridged_ring/PDB_rings/PHE2BCH_pairs/lig_menu/AaaaA_only4.csv"
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

    with open("phe2bch.err", "w") as log:
        log.write("Work starting \n")
    with open("phe2bch.log", "w") as log:
        log.write("Working starting \n")

    pdb_dir = "/home/chengyj/kinase_work/dataset/Bridged_ring/PDB_rings/PHE2BCH_pairs/pdb_dataset/pdb/"
    for key in tqdm(ligands_smi):
        if ligands_smi[key][1]:
            for pdbid in ligands_smi[key][2]:
                with open("phe2bch.log", "a") as log:
                    log.write(f"{pdbid}\n")
                pdb_file_path = pdb_dir + f"pdb{pdbid}.ent"

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

                try:
                    phe2bch_topdb(
                        ligands_smi[key][0], lig_Block, f"{key}_{pdbid}_bch.pdb"
                    )
                except Exception as ex:
                    with open("phe2bch.err", "a") as log:
                        log.write(f"Error of {key}_{pdbid}: {ex}" + "\n")
