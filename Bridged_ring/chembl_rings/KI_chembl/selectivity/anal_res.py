from pathlib import Path
from math import log
import numpy as np 
p=Path(".")
dat = sorted(list(p.glob("**/*selectivity.csv")))

all_data = {}
for item in dat:
    with item.open() as f:
        all_data[item.parts[-1].split(".")[0]] = np.array([log(float(item.split()[1])) for item in f.read().split("\n") ])
print()
for key in all_data:
    print(f"{key:>30} {all_data[key].size:>8} {np.mean(all_data[key]): >8.2f} {np.std(all_data[key]): >3.2f} {np.percentile(all_data[key],25): >3.2f} {np.percentile(all_data[key],50): >3.2f} {np.percentile(all_data[key],75): >3.2f}")
