from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, rdFMCS
from rdkit.Chem import rdRGroupDecomposition as rdRGD
from tqdm import tqdm

# take off warnings
RDLogger.DisableLog('rdApp.*')

def rename_atom(atom_to_rename,pdbinfo_ref):
    name = str(atom_to_rename.GetSymbol())+str(atom_to_rename.GetIdx()+1)
    if len(name) == 2:
        name = " " + name + " "
    elif len(name) == 3:
        name = " "+name
    elif len(name) == 4:
        pass
    pdbinfo = atom_to_rename.GetPDBResidueInfo()
    if pdbinfo:
        Chem.AtomMonomerInfo.SetName(pdbinfo,name)
        Chem.AtomPDBResidueInfo.SetAltLoc(pdbinfo," ")
        Chem.AtomPDBResidueInfo.SetResidueNumber(pdbinfo,1)
    else:
        Chem.Atom.SetMonomerInfo(atom_to_rename,pdbinfo_ref)
        pdbinfo = atom_to_rename.GetPDBResidueInfo()
        Chem.AtomMonomerInfo.SetName(pdbinfo,name)
    return name


def exctract_ligand_from_pdb(pdb_path,lig_id,smiles,outfile):
    # 抽取pdbfile中的ligand部分,获得rdkit可读字符串形式的pdbBlock
    # 并加氢重命名后以pdb格式储存在outfile中
    mol_from_smi = Chem.MolFromSmiles(smiles) 
    num = mol_from_smi.GetNumAtoms()
    ligand_lines = []
    with open(pdb_path) as ent:
        for line in ent:
            if line[0:6] == "HETATM" and line[17:20] == lig_id:
                ligand_lines.append(line)
    # print("".join(ligand_lines))
    ligand_start = 0
    ligand_end = len(ligand_lines)
    for i in range(1,len(ligand_lines)):
        if ligand_lines[i][17:26] != ligand_lines[i-1][17:26]:
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
    if ligand_end - ligand_start <= 2*num:
        lig_Block="".join(ligand_lines[ligand_start:ligand_end])
    else:
        lig_Block="".join(ligand_lines[ligand_start:ligand_start+num])
    
    # 加氢和重命名原子
    lig_pdb = Chem.MolFromPDBBlock(lig_Block)
    try:
        lig_pdb_refined = AllChem.AssignBondOrdersFromTemplate(mol_from_smi,lig_pdb)
    except:
        lig_pdb_refined = lig_pdb
    lig_pdb = Chem.AddHs(lig_pdb_refined,addCoords=True)
    atom_ref = lig_pdb.GetAtomWithIdx(1).GetPDBResidueInfo()
    for atom in lig_pdb.GetAtoms():
        rename_atom(atom,atom_ref)
    Chem.MolToPDBFile(lig_pdb,outfile)
    return lig_Block

def get_connect_atom_in_core(orimol,sidechain_atom_idxes):
    # mol is the whole molecule, the sidechain_atom_idx is the a tuple of the atom idx of a sidechain;
    # this function will get the atom connected with the sidechain in the core of the orimol.
    connect_atoms_in_core = []
    for idx in sidechain_atom_idxes:
        for atom in orimol.GetAtomWithIdx(idx).GetNeighbors():
            if atom.GetIdx() not in sidechain_atom_idxes:
                connect_atoms_in_core.append(atom.GetIdx())
    return connect_atoms_in_core

def get_connect_atom_in_sidechain(orimol,core):
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


def phe2bch_with_smiles(smi):
    # 用于将给定配体中的*间位苯环骨架*替换为*螺旋桨烷骨架*的脚本
    # 用smi指定ligand的结构

    # 读入分子
    mol = Chem.MolFromSmiles(smi)
    # 识别间位取代部分，* 匹配任意原子；'!:' 非芳香键； *!:ccc(!:*)ccc 含间位取代的苯环
    patt = Chem.MolFromSmarts("*!:c1cc(!:*)ccc1")
    match = mol.GetSubstructMatches(patt)
    
    # 识别两个间位取代基与苯环的连接原子,识别需要被替换苯环原子
    connect_atoms = get_connect_atom_in_sidechain(mol,match[0])
    core = set(match[0])-set(connect_atoms)

    # 记录第一部分的原子数目
    mol1_Atom_num = mol.GetNumAtoms()
    # 引入螺旋桨烷，并将其与原分子置于一处，使得其编号从·mol1_Atom_num+0·开始
    re_pat = Chem.MolFromSmiles("C12CCCC(C2)C1")
    mid_mol = Chem.CombineMols(mol,re_pat)
    # 开始成键
    mid_edited = Chem.EditableMol(mid_mol)
    # 螺旋桨烷烃的桥头碳分别为mol1_Atom_num和mol1_Atom_num+4
    mid_edited.AddBond(connect_atoms[0],mol1_Atom_num,order=Chem.rdchem.BondType.SINGLE)
    mid_edited.AddBond(connect_atoms[1],mol1_Atom_num+4,order=Chem.rdchem.BondType.SINGLE)
    # 删除原来的苯环
    # 每删除一个原子后原子序号都会发生改变，所以需要倒序删除。
    phenyl_atoms = list(core)
    phenyl_atoms.sort(reverse=True)
    for i in phenyl_atoms:
        mid_edited.RemoveAtom(i)
    BCHep_mol = mid_edited.GetMol()
    return BCHep_mol  

def getpdb(refmol,inpmol,pdbfile):
    
    bonded_conf = refmol.GetConformer()
    conf_res = inpmol.GetConformer()
    inpmol_prop = Chem.rdForceFieldHelpers.MMFFGetMoleculeProperties(inpmol)
    if inpmol_prop:
        ff_mcs = Chem.rdForceFieldHelpers.MMFFGetMoleculeForceField(inpmol,inpmol_prop)
    else:
        ff_mcs = Chem.rdForceFieldHelpers.UFFGetMoleculeForceField(inpmol)

    core = Chem.MolFromSmiles("C12CCCC(C2)C1")
    side_chains, unmatched = rdRGD.RGroupDecompose([core],[inpmol])
    for R in ['R1','R2']:
        sidechain = side_chains[0][R]
        sidechain_conf = sidechain.GetConformer()

        mcs_ref2mid = rdFMCS.FindMCS([refmol, sidechain],timeout=3,bondCompare=rdFMCS.BondCompare.CompareAny)
        for i, j in zip(refmol.GetSubstructMatch(mcs_ref2mid.queryMol),sidechain.GetSubstructMatch(mcs_ref2mid.queryMol)):
            sidechain_conf.SetAtomPosition(j, bonded_conf.GetAtomPosition(i))
        
        mcs_mid2inp = rdFMCS.FindMCS([inpmol, sidechain],timeout=3,bondCompare=rdFMCS.BondCompare.CompareOrderExact)
        for i, j in zip(sidechain.GetSubstructMatch(mcs_mid2inp.queryMol),inpmol.GetSubstructMatch(mcs_mid2inp.queryMol)):
            ff_mcs.AddFixedPoint(j)
            conf_res.SetAtomPosition(j, sidechain_conf.GetAtomPosition(i))
        
        mcs1 = rdFMCS.FindMCS([refmol, sidechain],timeout=3,bondCompare=rdFMCS.BondCompare.CompareAny)
        mcs3 = rdFMCS.FindMCS([inpmol, sidechain],timeout=3,bondCompare=rdFMCS.BondCompare.CompareOrderExact)
        sidechain_ref = refmol.GetSubstructMatch(mcs1.queryMol)
        sidechain_gen = inpmol.GetSubstructMatch(mcs3.queryMol)
        atom_connect_sidechain_core_ref = get_connect_atom_in_core(refmol, sidechain_ref)[0]
        atom_connect_sidechain_core_gen = get_connect_atom_in_core(inpmol, sidechain_gen)[0]
        ff_mcs.AddFixedPoint(atom_connect_sidechain_core_gen)
        conf_res.SetAtomPosition(atom_connect_sidechain_core_gen,bonded_conf.GetAtomPosition(atom_connect_sidechain_core_ref))
     
    for i in range(10):
        try:
            ff_mcs.Minimize()
        except:
            pass
    inpmol = Chem.AddHs(inpmol,addCoords=True)
    Chem.MolToPDBFile(inpmol, pdbfile)
    return inpmol

def phe2bch_topdb(smi0,refpdb,name):
    mol = phe2bch_with_smiles(smi0)
    smi = Chem.MolToSmiles(mol)
    mol = Chem.MolFromSmiles(smi)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol,maxAttempts=5000)
    mol = Chem.RemoveAllHs(mol)

    structure_from_pdb = Chem.MolFromPDBBlock(refpdb)
    try:
        structure_from_pdb.GetConformer()
    except:
        structure_from_pdb = Chem.MolFromPDBBlock(refpdb,sanitize=False,flavor=1)
    
    try:
        structure_from_pdb = Chem.RemoveAllHs(structure_from_pdb)
    except:
        structure_from_pdb = Chem.RemoveHs(structure_from_pdb,implicitOnly=True)
    # some ligand can not remove all H
    
    mol_tem = Chem.MolFromSmiles(smi0)
    try:
        structure_refine = AllChem.AssignBondOrdersFromTemplate(mol_tem,structure_from_pdb)
    except:
        structure_refine = structure_from_pdb
    getpdb(structure_refine,mol,f"{name}")

if __name__ == "__main__":
    ligands_smi = {}
    with open("/home/chengyj/kinase_work/dataset/Bridged_ring/PDB_rings/PHE2BCH_pairs/lig_menu/AaaaA_only4.csv") as AaaaA:
        for line in AaaaA:
            info = line.split()
            ligands_smi[info[1]] = [info[0]]
    # print(1)
    with open("/home/chengyj/kinase_work/dataset/Bridged_ring/PDB_rings/PHE2BCH_pairs/pdb_dataset/cc-to-pdb.tdd") as ligs:
        for line in ligs:
            info = line.split()
            if info[0] in ligands_smi.keys():
                ligands_smi[info[0]].append(info[1:])
    # print(2)

    with open("phe2bch.log","w") as log:
        log.write("Work starting \n")
    for key in tqdm(ligands_smi):
        # print(key)
        try:
            for pdbid in ligands_smi[key][1]:
                pdb_file_path=f"/home/chengyj/kinase_work/dataset/Bridged_ring/PDB_rings/PHE2BCH_pairs/pdb_dataset/pdb/pdb{pdbid}.ent"
                lig_Block = exctract_ligand_from_pdb(pdb_file_path,key,ligands_smi[key][0],f"{key}_{pdbid}_phe.pdb")
                try:
                    phe2bch_topdb(ligands_smi[key][0],lig_Block,f"{key}_{pdbid}_bch.pdb")
                except Exception as ex:
                    with open("phe2bch.log","a") as log:
                        log.write(f"Error of {key}_{pdbid}: {ex}"+"\n")
        except:
            pass
