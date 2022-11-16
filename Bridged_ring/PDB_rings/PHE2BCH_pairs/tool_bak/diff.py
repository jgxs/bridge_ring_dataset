import sys 
import re
class atom:
    def __init__(self,line):
        info = line.split()
        self.line = line
        self.name = info[2]
        self.x = float(line[31:39])
        self.y = float(line[39:48])
        self.z = float(line[48:56])
        self.cor=(self.x,self.y,self.z)
        self.type = re.sub(r'[0-9]',"",self.name)
    def getdist(self,atom):
        return ((atom.x-self.x)**2+(atom.y-self.y)**2+(atom.z-self.z)**2)**0.5 

def get_diff_atoms(lig1,lig2):
    diff_atoms_real = []
    diff_atoms_change = []
    for item in lig1:
        isdiff_s = True
        isdiff_c = True
        if item in lig2 and lig1[item].type == lig2[item].type:
            isdiff_s = False
            isdiff_c =False
        else:
            for cor in lig2:
                if lig1[item].type==lig2[cor].type:
                    if lig1[item].getdist(lig2[cor]) < 0.5:
                        isdiff_s = False
                    if lig1[item].getdist(lig2[cor]) < 0.01:
                        isdiff_c = False
        if isdiff_s:
            diff_atoms_real.append(lig1[item].name)
        if isdiff_c and not isdiff_s:
            diff_atoms_change.append(lig1[item].line)
    return [diff_atoms_real,diff_atoms_change]

lig=[{}]
with open(sys.argv[1]) as ligspdb:
    for line in ligspdb:
        if line[0:6] == "HETATM": 
            lig[-1][atom(line).cor] = atom(line)
            #print(atom(line).cor)
        else:
            lig.append({})



diff=[[],[],[],[]]
diff[0],diff[2] = get_diff_atoms(lig[0],lig[1])
diff[1],diff[3] = get_diff_atoms(lig[1],lig[0])
with open("tmp.pdb","w") as out:
    for key in lig[1]:
        if lig[1][key].line not in diff[3]:
            out.write(lig[1][key].line)
        else:
            line = lig[1][key].line
            dist = 10000.0
            for cor in lig[0]:
                if lig[0][cor].getdist(lig[1][key]) < dist:
                    dist = lig[0][cor].getdist(lig[1][key])
                    line = line[0:27]+lig[0][cor].line[27:]
            out.write(line)


print(str(diff[0])[2:-2].replace("', '",",")+"_"+str(diff[1])[2:-2].replace("'",""))
#print(str(diff[2])+"_"+str(diff[3]))
# ['C4', 'C5', 'C6', 'H2', 'H3', 'H4', 'H8', 'H9', 'H12', 'H13', 'H14', 'H37', 'H38', 'H39']_['C34', 'C35', 'C36', 'H7', 'H8', 'H9', 'H25', 'H26', 'H27', 'H42', 'H47', 'H48', 'H49', 'H50']]