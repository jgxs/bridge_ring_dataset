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

# when sidechain is a single C
```
def get_correct_match(mol_whole, mcs, arranged_atoms, core_mol, debug):
    match_to_use = False
    correct_matches = mol_whole.GetSubstructMatches(mcs.queryMol)
    print(f"test len{correct_matches}")
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
        core_head, side_head, core_all = get_connecter(mol_whole,core_mol)
        for match_correct in correct_matches:
            if len(match_correct) == 1:
                # print(set(match_correct))
                # print(set(match_correct) & side_head)
                print(set(match_correct) & core_all)
                print(set(match_correct) & core_all != set())
            if set(match_correct) & core_head == set() and set(match_correct) & arranged_atoms == set() and set(match_correct) & core_all == set():
                if len(match_correct) == 1:
                    print("hello world")
                match_to_use = match_correct
    else:
        if debug:
            print(f"begin4 {correct_matches}")
        match_to_use = correct_matches[0]
    if debug:
        pass
    
    #    print(f"end {match_to_use}")
    if match_to_use:
        return match_to_use
    else:
        print("i really do not want")
        return correct_matches[0]
```
