from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors
from scipy.stats import norm
import matplotlib.pyplot as plt
import numpy as np 

from multiprocessing import Pool, Lock

from pathlib import Path
import argparse


def calFsp3(line_in_smis):
    smi = line_in_smis.split()[1]
    try:
        fsp3 = rdMolDescriptors.CalcFractionCSP3(Chem.MolFromSmiles(smi))
        print(f" {fsp3} ", end=" ")
    except:
        pass
    return True
"""
def init(l):
    global lock
    lock = l
"""

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Plot TI figure")
    parser.add_argument(
        "--data",
        dest='data',
    )
    smis = Path(parser.parse_args().data)
    with smis.open() as info:
        allmols = info.read().split("\n")[0:-1]
    
    lock = Lock()
    with Pool(4) as p:
        p.map(calFsp3, allmols)
    """
    flat = allFsp
    mu = np.mean(flat)
    print(max(flat))
    print(min(flat))
    sigma = np.std(flat)
    num_bins = 19
    n, bins, patches = plt.hist(flat, num_bins, density=True,stacked=True,facecolor='blue', alpha=0.5)
    print(bins)
    y = norm.pdf(bins,mu,sigma)
    plt.plot(bins, y, 'r--') #绘制y的曲线
    plt.text(0.1,2.5,f"$\mu = {mu:.2f},\ \sigma = {sigma:.2f}$")
    plt.xlabel('Fsp3 of compounds') #绘制x轴
    plt.ylabel('Probability') #绘制y轴
    plt.title(f'Distribution of Fsps')#中文标题 u'xxx' 
    plt.subplots_adjust(left=0.15)#左边距 
    plt.savefig(f"{output}.png")
    plt.close()
    plt.boxplot([flat],sym="")
    plt.savefig("test.png")
    """
