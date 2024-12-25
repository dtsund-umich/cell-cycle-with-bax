import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import scipy.stats
import os
import sys

inkdata = []
baxdata = []

ignore_steady = False
if len(sys.argv) >= 2 and sys.argv[1] == "-i":
    ignore_steady = True

for i in range(2500):
    if ignore_steady and os.path.exists(str(i)+"/STEADY.txt"):
        continue
    
    inktext = open(str(i)+"/Ink4.txt",'r').readlines()
    baxtext = open(str(i)+"/BAX.txt",'r').readlines()
    
    inktemp = []
    baxtemp = []
    
    for line in inktext:
        inktemp.append(float(line.strip().split()[1]))
    for line in baxtext:
        baxtemp.append(float(line.strip().split()[1]))
    
    inkdata.append(max(inktemp[len(inktemp)*4//5:-1]))
    baxdata.append(max(baxtemp[len(baxtemp)*4//5:-1]))
    #inkdata.append(sum(inktemp[len(inktemp)*4//5:-1])/(len(inktemp)*4//5))
    #baxdata.append(sum(baxtemp[len(baxtemp)*4//5:-1])/(len(baxtemp)*4//5))

plt.xlabel("Peak Ink4")
plt.ylabel("Peak BAX")

print(scipy.stats.pearsonr(inkdata, baxdata))
plt.scatter(inkdata,baxdata)
plt.show()
