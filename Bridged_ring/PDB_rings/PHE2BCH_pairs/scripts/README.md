# Some tips about the usage of rdkit

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



