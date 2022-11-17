# Bioisostere Dataset of PHE-BCH molecular pair

Here is my work of the Bioisostere Dataset of PHE-BCH molecular pair.

## Inspiration

Inspirated by the [Synthesis of Meta-Substituted Arene Bioisosteres from [3.1.1]Propellane.](https://doi.org/10.1038/s41586-022-05290-z), we want to prove the potient of the PHE-BCH molecular pair.

## Workflow

First, use `lig_AaaaA_sel.py` to select the ligands with meta-substituted arene group. In order to the group replacement in later step, we only choose the phenyl with two substitued groups. There are 942 ligands in PDB as required.

Then, use `Phe2BCH.py` to generate the complex structure of the ligand with BCHep. We constuct 1177 complexes of 864 ligands. Here, we add the H atoms of the ligands through rdkit with the help of the smiles.

Then, we plan to use plop to rescore the complex and calculate the interaction energy of the difference molecular pairs.