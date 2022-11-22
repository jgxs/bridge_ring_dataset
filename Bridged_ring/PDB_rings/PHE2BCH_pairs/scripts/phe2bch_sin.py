from Phe2BCH import *


ligands_smi = {}
with open(
    "/home/chengyj/kinase_work/dataset/Bridged_ring/PDB_rings/PHE2BCH_pairs/lig_menu/AaaaA_only3.csv"
) as AaaaA:
    for line in AaaaA:
        info = line.split()
        ligands_smi[info[1]] = [info[0]]
# inp = "0FS_3vc4" segmentation fault

AaaaA_core = Chem.MolFromSmarts("*!:c1cc(!:*)ccc1")
bch_core = Chem.MolFromSmiles("C12CCCC(C2)C1")
def get_connecter(mol,patt):
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

def get_correct_match(mol_whole, mcs, arranged_atoms, core_mol, debug):
    correct_matches = mol_whole.GetSubstructMatches(mcs.queryMol)
    if debug:
        print(f"begin {correct_matches}")
        print(arranged_atoms)
    if not arranged_atoms and len(correct_matches) > 1:
        if debug:
            print(f"begin2 {correct_matches}")
        # arranged_atoms is not empty means it is the first bigger group to match
        # however, the len(correct_matches) > 1 means the bigger group is not 
        # big enough to match only one part of the refmol
        # In that situation, we need to select the correct part of the molecule.
        core_head, side_head = get_connecter(mol_whole,core_mol)
        print(core_head)
        print(side_head)
        for match_correct in correct_matches:
            print(set(match_correct) & core_head)
            print(set(match_correct) & side_head)
            if set(match_correct) & core_head == set() and set(match_correct) & side_head != side_head:
                match_to_use = match_correct
    elif len(correct_matches) > 1 and arranged_atoms:
        if debug:
            print(f"begin3 {correct_matches}")
        for match_correct in correct_matches:
            if set(match_correct) & arranged_atoms != set():
                pass
            else:
                match_to_use = match_correct
                break 
    else:
        if debug:
            print(f"begin4 {correct_matches}")
        match_to_use = correct_matches[0]
    if debug:
        print(f"end {match_to_use}")
    return match_to_use

def get_connecter(mol,patt):
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

def getpdb_debug(refmol, inpmol, pdbfile):

    bonded_conf = refmol.GetConformer()
    conf_res = inpmol.GetConformer()
    inpmol_prop = Chem.rdForceFieldHelpers.MMFFGetMoleculeProperties(inpmol)
    if inpmol_prop:
        ff_mcs = Chem.rdForceFieldHelpers.MMFFGetMoleculeForceField(inpmol, inpmol_prop)
    else:
        ff_mcs = Chem.rdForceFieldHelpers.UFFGetMoleculeForceField(inpmol)

    core = Chem.MolFromSmiles("C12CCCC(C2)C1")
    side_chains, unmatched = rdRGD.RGroupDecompose([core], [inpmol])
    if side_chains[0]["R1"].GetNumAtoms()>side_chains[0]["R2"].GetNumAtoms():
        side_chains_reorder = [side_chains[0]["R1"],side_chains[0]["R2"]]  
    else:
        side_chains_reorder = [side_chains[0]["R2"],side_chains[0]["R1"]]
    
    ref_arranged_atoms = set()
    inp_arranged_atoms = set()
    for sidechain in side_chains_reorder:
        # Chem.SanitizeMol(sidechain)
        # print(ref_arranged_atoms)
        # print(inp_arranged_atoms)
        for item in sidechain.GetAtoms():
            print(item.GetSymbol())
        sidechain_conf = sidechain.GetConformer()
        corr = {}
        mcs1 = rdFMCS.FindMCS(
            [refmol, sidechain], timeout=3, bondCompare=rdFMCS.BondCompare.CompareAny
        )
        mcs3 = rdFMCS.FindMCS(
            [inpmol, sidechain],
            timeout=3,
            bondCompare=rdFMCS.BondCompare.CompareOrderExact,
        )

        sidechain_ref = get_correct_match(refmol,mcs1,ref_arranged_atoms,AaaaA_core,False)
        sidechain_gen = get_correct_match(inpmol,mcs3,inp_arranged_atoms,bch_core,False)
        
        atom_connect_sidechain_core_ref = get_connect_atom_in_core(refmol, sidechain_ref)[0]
        atom_connect_sidechain_core_gen = get_connect_atom_in_core(inpmol, sidechain_gen)[0]

        ff_mcs.AddFixedPoint(atom_connect_sidechain_core_gen)
        conf_res.SetAtomPosition(atom_connect_sidechain_core_gen,bonded_conf.GetAtomPosition(atom_connect_sidechain_core_ref),)
        
        mcs_ref2mid = rdFMCS.FindMCS([refmol, sidechain], timeout=3, bondCompare=rdFMCS.BondCompare.CompareAny)
        match_to_correct_ref = get_correct_match(refmol,mcs_ref2mid,ref_arranged_atoms,AaaaA_core,False)
        print(f"the result: {sidechain.GetSubstructMatches(mcs_ref2mid.queryMol)}")
        for i, j in zip(sidechain.GetSubstructMatch(mcs_ref2mid.queryMol),match_to_correct_ref,):
            sidechain_conf.SetAtomPosition(i, bonded_conf.GetAtomPosition(j))
            print(f"set atom idx {i}")
            ref_arranged_atoms.add(j)
        # print(ref_arranged_atoms)
        mcs_mid2inp = rdFMCS.FindMCS(
            [inpmol, sidechain],
            timeout=3,
            bondCompare=rdFMCS.BondCompare.CompareOrderExact,
        )
        match_to_correct_inp = get_correct_match(inpmol,mcs_mid2inp,inp_arranged_atoms,bch_core,True)
        print(match_to_correct_inp)
        print(f"sidechain {sidechain.GetSubstructMatches(mcs_mid2inp.queryMol)}")
        for i, j in zip(
            sidechain.GetSubstructMatch(mcs_mid2inp.queryMol),
            match_to_correct_inp
        ):  
            print(f"{i}      i")
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

def phe2bch_topdb_debug(smi0, refpdb, name):
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
    return getpdb_debug(structure_refine, mol, f"{name}")

# print("test")
inp = "0OO_4eh7"
# inp = "0RZ_4f1q"
inp = "07Q_3tv7"
inp = "B74_6xic"
key = inp.split("_")[0]
pdb_id = inp.split("_")[1]
lig_smi = ligands_smi[key][0]
# print(lig_smi)
pdb_path = f"/home/chengyj/kinase_work/dataset/Bridged_ring/PDB_rings/PHE2BCH_pairs/pdb_dataset/pdb/pdb{pdb_id}.ent"
print(pdb_path)
lig_Block = exctract_ligand_from_pdb(pdb_path, key, lig_smi, f"{key}_{pdb_id}_sin.pdb")
test = phe2bch_topdb_debug(lig_smi, lig_Block, f"{key}_{pdb_id}_bch.pdb")
