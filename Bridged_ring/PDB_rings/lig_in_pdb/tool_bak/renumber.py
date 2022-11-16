import argparse
import re

parser = argparse.ArgumentParser(description="ref_lib")
parser.add_argument(
    "--lib",
    dest='lib',
    )
parser.add_argument(
    "--pdb",
    dest='pdb'
)
lib = parser.parse_args().lib
pdb = parser.parse_args().pdb

with open(lib) as info:
    allinfo = info.read()
ref_atom_name = [item.split()[0][1:-1] for item in allinfo.split("!entry.")[1].split("\n")[1:-1] ]

i=0
with open(pdb) as pdbfile:
    for line in pdbfile:
        element = re.sub(r'[0-9]',"",ref_atom_name[i])
        if len(element) < 2:
            print(f"{line[0:12]} {ref_atom_name[i]:<3} {line[17:-1]}")
        else:
            print(f"{line[0:12]}{ref_atom_name[i]:<4} {line[17:-1]}")
        i+=1
        if i == len(ref_atom_name):
            break
 