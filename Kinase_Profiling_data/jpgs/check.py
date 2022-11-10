from rdkit import Chem
from rdkit.Chem import Draw

smis = []
with open("ligands.smi") as ligands:
    for line in ligands:
        smis.append([line.split()[0],line.split()[-1]])

for item in smis:
    m = Chem.MolFromSmiles(item[0])
    print(item[1].split('.')[0])
    #Draw.MolToImage(m,size=(1400,1400))
    Draw.MolToFile(m, f"jpg_from_smi/{item[1].split('.')[0]}.png", size=(1500, 1500),bondLineWidth=5.0)
