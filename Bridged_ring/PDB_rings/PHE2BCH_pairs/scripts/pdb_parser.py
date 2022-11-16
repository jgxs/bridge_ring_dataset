from pathlib import Path

class Atom:
    # parse a line in pdb file and generate a Atom contain information of it.
    def __init__(self,atomline):
        atominfo = atomline.split()
        self.name = atominfo[2]
        self.element = atominfo[2][0]
        self.xyz = [float(atomline[30:38]),float(atomline[38:46]),float(atomline[46:54])]

    def get_dist(self,atom):
        return sum([(self.xyz[i]-atom.xyz[i])**2 for i in range(3)])**0.5

class Ligand:
    def __init__(self,pdbpath):
        self.lines = []
        self.atoms={}
        i = 0
        with pdbpath.open() as pdbfile:
            for line in pdbfile:
                if line[0:6] == "HETATM":
                    self.lines.append(line)
                    self.atoms[i] = Atom(line)
                    i += 1
        self.chain_id = self.lines[0][21]
        self.number = i

class chain:
    def __init__(self,pdbfile,chain_id,target_path) -> None:
        pdbpath = Path(pdbfile)
        self.pdbid = pdbpath.name.replace(".pdb","")
        self.PDBBlock = []
        with pdbpath.open() as pdblines:
            for line in pdblines:
                if line[0:4] == "ATOM" and line[21] == chain_id:
                    self.PDBBlock.append(line)

        with open(f"{target_path}/{self.pdbid}_sin.pdb","w") as output:
            for item in self.PDBBlock:
                output.write(item)
            output.write("TER\n")


if __name__ == "__main__":
    Ligand(Path("/home/chengyj/kinase_work/dataset/Bridged_ring/PDB_rings/lig_in_pdb/lig_bch_complex/62Y/4cxu/62Y_4cxu_phe.pdb"))


