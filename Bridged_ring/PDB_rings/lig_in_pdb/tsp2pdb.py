from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS

def getpdb(refmol,inpmol,pdbfile):
    mcs = rdFMCS.FindMCS([refmol, inpmol],timeout=3)
    bonded_conf = refmol.GetConformer()
    conf_res = inpmol.GetConformer()
    inpmol_prop = Chem.rdForceFieldHelpers.MMFFGetMoleculeProperties(inpmol)
    ff_mcs = Chem.rdForceFieldHelpers.MMFFGetMoleculeForceField(inpmol,inpmol_prop)
    for i, j in zip(refmol.GetSubstructMatch(mcs.queryMol), inpmol.GetSubstructMatch(mcs.queryMol)):
        ff_mcs.AddFixedPoint(j)
        conf_res.SetAtomPosition(j, bonded_conf.GetAtomPosition(i))
    try:
        ff_mcs.Minimize()
    except:
        pass
    Chem.MolToPDBFile(inpmol, pdbfile)
    return inpmol
#print(rdBase.rdkitVersion)
with open("inputs") as infor:
    inf = infor.read().split()

smi1 = inf[2]
smi2 = inf[3]
pdb = "lig_template.pdb"
out = inf[1]
print(smi1)
structure_from_pdb = Chem.MolFromPDBFile(pdb)
mol_tem = Chem.MolFromSmiles(smi1)
structure_refine = AllChem.AssignBondOrdersFromTemplate(mol_tem,structure_from_pdb)
names= ["N","C"]
smis = [smi1, smi2]
for item in range(2):
    mol = Chem.MolFromSmiles(smis[item])
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    mol = Chem.RemoveAllHs(mol)
    getpdb(structure_refine,mol,f"{out}_{names[item]}.pdb")
