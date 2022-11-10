from scipy.stats import norm,mode
import matplotlib.pyplot as plt
import numpy as np 
from pathlib import Path
p=Path(".")
dat = sorted(list(p.glob("**/*fsp3")))
all_data = {}
for item in dat:
    with item.open() as f:
        all_data[item.parts[-1].split(".")[0]] = np.array([float(item) for item in f.read().split() ])
print()
for key in all_data:
    print(f"{key:>10} {all_data[key].size:>8} {np.mean(all_data[key]):.3f} {np.std(all_data[key]):.3f} {np.median(all_data[key]):.3f}")

plt.figure(figsize(10,5))
plt.title('The distribution of F-sp3 of compounds in different database')
labels = all_data.keys()
plt.boxplot(list(all_values()),labels=labels,)