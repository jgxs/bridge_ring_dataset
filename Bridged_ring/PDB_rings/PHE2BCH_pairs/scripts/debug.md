# can not understand
```python
    mol = phe2bch_with_smiles(smi0)
    ## Totally can not understand
    # here must mol to smi to mol
    # otherwise, the Allchem.EmbedMolecule may be
    # Segmentation fault
    # the bug detail is the "0FS_3vc4"
    smi = Chem.MolToSmiles(mol)
    mol = Chem.MolFromSmiles(smi)
    mol = Chem.AddHs(mol)
```