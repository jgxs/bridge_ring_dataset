# Some tips about the usage of rdkit

## MolFromPDBBlock and MolFromPDBFile

### flavor

Sometimes, `MolFromPDBBlock` and `MolFromPDBFile` cannot get Mol object. A possible way to solve this problem is set the optional `flavor` as 1. 

The reason of this problem may be `PDBParser` will ignore some atoms when flavor are set to 0, the default value.
```C++
    // ptr is the atom line in pdb file 
  if ((flavor & 1) == 0) {
    // Ignore alternate locations of atoms.
    if (len >= 17 && ptr[16] != ' ' && ptr[16] != 'A' && ptr[16] != '1') {
      return;
    }
    // Ignore XPLOR pseudo atoms
    if (len >= 54 && !memcmp(ptr + 30, "9999.0009999.0009999.000", 24)) {
      return;
    }
    // Ignore NMR pseudo atoms
    if (ptr[12] == ' ' && ptr[13] == 'Q') {
      return;
    }
    // Ignore PDB dummy residues
    if (len >= 20 && !memcmp(ptr + 18, "DUM", 3)) {
      return;
    }
  }

```

### proximityBonding

When there is an 1,2,3-triazole group in the mol, the `MolFromPDBBlock` and `MolFromPDBFile` may connect the atoms uncorrectly. You can see the [bug in detail](https://github.com/rdkit/rdkit/issues/5778)

A simple way to avoid this problem is to set the proximityBonding = False, like

```python 
triazole_from_pdb = Chem.MolFromPDBFile("triazole.pdb",proximityBonding=False)
```
However, that only can handle the mol with few atoms. 

Here, display another way the aovid that problem: remove the wrong bond by comparing the length of bonds

```python
def check_triazole_pdb(mol_pdb,smi):
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
        for atom in mol_pdb.GetAtoms():
            atom_symbol = atom.GetSymbol()
            neighors = atom.GetNeighbors()
            if atom_symbol == "N" and len(neighors) > 3:
                N_idx = atom.GetIdx()
                dist_Idx = [(Point3D.Distance(conf.GetAtomPosition(item.GetIdx()),conf.GetAtomPosition(N_idx)), item.GetIdx()) for item in neighors]
        dist_Idx.sort()
        mol_edited.RemoveBond(N_idx,dist_Idx[-1][1])
        return mol_edited.GetMol()

```

## AtomPDBResidueInfo and AtomMonomerInfo

When load a Mol from a pdb file by `Chem.MolFromPDBFile()` or `Chem.MolFromPDBBlock()`, the information of residue will be stored in the [AtomPDBResidueInfo](https://www.rdkit.org/docs/source/rdkit.Chem.rdchem.html#rdkit.Chem.rdchem.AtomPDBResidueInfo) class. We can use scripts below to get and change the information.
```python
Info = atom.GetPDBResidueInfo()
Chem.AtomPDBResidueInfo.GetAltLoc(pdbinfo)
Chem.AtomPDBResidueInfo.SetResidueNumber(pdbinfo,1)
```
As for the H atoms added by the `Chem.AddHs()`, the return value of `Hs.GetPDBResidueInfo()` is `None`. One possible method to set these information is shown below. We can copy another atom's `AtomPDBResidueInfo` to it and change the value after it.
```python
Atom_noH_info = Mol_from_pdb.GetAtomWithIdx(1).GetPDBResidueInfo()
H_info = Chem.Atom.SetMonomerInfo(H_atom,Atom_noH_info)
pdbinfo = H_info.GetPDBResidueInfo()
Chem.AtomMonomerInfo.SetName(H_info,name)
```

## RemoveHs() and RemoveAllHs()

`RemoveHs()` can remove only the implicit Hs of the mol, which may be more useful than `RemoveAllHs()`

## AllChem.AssignBondOrdersFromTemplate(refmol, inpmol)

Template molecule should have no explicit hydrogens. A possible way to solve that problem is:
```python 
remol = Chem.RemoveAllHs(mol_from_smi)
lig_pdb_refined = AllChem.AssignBondOrdersFromTemplate(remol,lig_pdb)
```

