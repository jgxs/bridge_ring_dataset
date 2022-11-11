from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
from rdkit.Chem import rdRGroupDecomposition as rdRGD

def pdb_2_lig_block(pdb_path,lig_id):
    # 抽取pdbfile中的ligand部分,获得rdkit可读字符串形式的pdbBlock
    ligand_lines = []
    with open(pdb_path) as ent:
        for line in ent:
            if line[0:6] == "HETATM" and line[17:20] == lig_id:
                ligand_lines.append(line)
    ligand_num = len(ligand_lines)
    for i in range(1,len(ligand_lines)):
        if ligand_lines[i][17:26] != ligand_lines[i-1][17:26]:
            ligand_num = i
            break
        else:
            pass
    lig_Block="".join(ligand_lines[0:ligand_num])
    return lig_Block

def get_connect_atom_in_core(orimol,sidechain_atom_idxes):
    # mol is the whole molecule, the sidechain_atom_idx is the a tuple of the atom idx of a sidechain;
    # this function will get the atom connected with the .
    connect_atoms_in_core = []
    for idx in sidechain_atom_idxes:
        for atom in orimol.GetAtomWithIdx(idx).GetNeighbors():
            if atom.GetIdx() not in sidechain_atom_idxes:
                connect_atoms_in_core.append(atom.GetIdx())
    return connect_atoms_in_core

def get_connect_atom_in_sidechain(orimol,core):
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
    patt = Chem.MolFromSmarts("*!:ccc(!:*)ccc")
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
    ff_mcs = Chem.rdForceFieldHelpers.MMFFGetMoleculeForceField(inpmol,inpmol_prop)

    core = Chem.MolFromSmiles("C12CCCC(C2)C1")
    res, unmatched = rdRGD.RGroupDecompose([core],[inpmol])
    for R in ['R1','R2']:
        mcs1 = rdFMCS.FindMCS([refmol, res[0][R]],timeout=3,bondCompare=rdFMCS.BondCompare.CompareOrder)
        mcs3 = rdFMCS.FindMCS([inpmol, res[0][R]],timeout=3,bondCompare=rdFMCS.BondCompare.CompareOrderExact)
        sidechain_ref = refmol.GetSubstructMatch(mcs1.queryMol)
        sidechain_gen = inpmol.GetSubstructMatch(mcs3.queryMol)
        for i, j in zip(sidechain_ref,sidechain_gen):
            ff_mcs.AddFixedPoint(j)
            conf_res.SetAtomPosition(j, bonded_conf.GetAtomPosition(i))
        atom_connect_sidechain_core_ref = get_connect_atom_in_core(refmol, sidechain_ref)[0]
        atom_connect_sidechain_core_gen = get_connect_atom_in_core(inpmol, sidechain_gen)[0]
        ff_mcs.AddFixedPoint(atom_connect_sidechain_core_gen)
        conf_res.SetAtomPosition(atom_connect_sidechain_core_gen,bonded_conf.GetAtomPosition(atom_connect_sidechain_core_ref))
    for i in range(6):
        try:
            ff_mcs.Minimize()
        except:
            pass
    Chem.MolToPDBFile(inpmol, pdbfile)
    return inpmol

def phe2bch_topdb(smi0,refpdb,name):
    mol = phe2bch_with_smiles(smi0)
    smi = Chem.MolToSmiles(mol)
    mol = Chem.MolFromSmiles(smi)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    mol = Chem.RemoveAllHs(mol)

    structure_from_pdb = Chem.MolFromPDBBlock(refpdb)
    structure_from_pdb = Chem.RemoveAllHs(structure_from_pdb)
    mol_tem = Chem.MolFromSmiles(smi0)
    structure_refine = AllChem.AssignBondOrdersFromTemplate(mol_tem,structure_from_pdb)
    getpdb(structure_refine,mol,f"{name}")

if __name__ == "__main__":
    ligands_smi = {}
    with open("/home/chengyj/kinase_work/dataset/Bridged_ring/PDB_rings/lig_in_pdb/lig_AaaaA_test3.csv") as AaaaA:
        for line in AaaaA:
            info = line.split()
            ligands_smi[info[1]] = [info[0]]
    print(1)
    with open("/home/chengyj/kinase_work/dataset/Bridged_ring/PDB_rings/lig_in_pdb/pdb_dataset/cc-to-pdb.tdd") as ligs:
        for line in ligs:
            info = line.split()
            if info[0] in ligands_smi.keys():
                ligands_smi[info[0]].append(info[1:])
    print(2)
    for key in ligands_smi:
        print(key)
        try:
            for pdbid in ligands_smi[key][1]:
                print(pdbid)
                """
                lig_lines = []                
                with open(f"/home/chengyj/kinase_work/dataset/Bridged_ring/PDB_rings/lig_in_pdb/pdb_dataset/pdb/pdb{pdbid}.ent") as pdbfile:
                    for line in pdbfile:
                        if line[0:6] == "HETATM" and line[17:20] == key:
                            lig_lines.append(line)
                ligand_num = len(lig_lines)
                for i in range(1,len(lig_lines)):
                    if lig_lines[i][17:26] != lig_lines[i-1][17:26]:
                        ligand_num = i
                        break
                    else:
                        pass
                lig_Block = "".join(lig_lines[0:ligand_num])
                """
                pdb_file_path=f"/home/chengyj/kinase_work/dataset/Bridged_ring/PDB_rings/lig_in_pdb/pdb_dataset/pdb/pdb{pdbid}.ent"
                lig_Block = pdb_2_lig_block(pdb_file_path,key)
                with open(f"{key}_{pdbid}_phe.pdb",'w') as phe:
                    for item in lig_Block:
                        phe.write(item)
                phe2bch_topdb(ligands_smi[key][0],lig_Block,f"{key}_{pdbid}_bch.pdb")
        except:
            pass
