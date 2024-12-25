import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import os
import sys


#380, 381, ... 399
#...
#20, 21, ... 39
#0, 1, ..., 19

#^^^this, except 50x50 instead of 20x20
nums = []

for i in reversed(range(0,50)):
    nums += range(50*i,50*i+50)

xlabel = ""
ylabel = ""
xticks = []
yticks = []

if len(sys.argv) > 1:
    infile = open(sys.argv[1],"r").readlines()
    xlabel = infile[0].split()[0].strip()
    ylabel = infile[1].split()[0].strip()
    xlower = float(infile[0].split()[1].strip())
    xupper = float(infile[0].split()[2].strip())
    ylower = float(infile[1].split()[1].strip())
    yupper = float(infile[1].split()[2].strip())
    xticks = ["",str(xlower), str((4*xlower+xupper)/5), str((3*xlower+2*xupper)/5), str((2*xlower+3*xupper)/5), str((xlower+4*xupper)/5),str(xupper)]
    yticks = ["",str(yupper), str((4*yupper+ylower)/5), str((3*yupper+2*ylower)/5), str((2*yupper+3*ylower)/5), str((yupper+4*ylower)/5),str(ylower)]
else:
    xlabel = "E6"
    ylabel = "E7"
    xticks = ["","0","2","4","6","8","10"]
    yticks = ["","3","2.4","1.8","1.2","0.6","0"]
    
"""
mask = False
maskdata = []
if len(sys.argv) > 2:
    mask = True
    infile = open(sys.argv[2],"r").readlines()
    for i in infile:
        maskdata.append(i.strip())
"""
    

#Get all the BAX peaks
peaks = []
maskdata = []
for i in nums:
    data = []

    if os.path.exists(str(i)+"/STEADY.txt"):
        maskdata.append(1)
        peaks.append(np.nan)
    else:
        maskdata.append(0)
        filedata = open(str(i)+"/BAX.txt",'r').readlines()
        for line in filedata:
            data.append(float(line.strip().split()[1]))
        #If the data starts at t=0, cut out initial transient behavior
        if filedata[0].split()[0] == "0.0":
            data = data[int(len(data)*4/5):-1]
        peak = max(data)
        #peak = float(open(str(i)+"/BAX.txt",'r').readlines()[0].strip().split()[1])
        peaks.append(peak)
        
        """
        if peak < 20:
            peaks.append(peak)
        else:
            peaks.append(20)
        """
    
    #lines = open(str(i)+"/Me.txt",'r').readlines()
    #temp = []
    #for line in lines:
    #    temp.append(float(line.split()[1].strip()))
    #
    #if max(temp[len(temp)/2:]) < 2.5:
    #    peaks.append(max(temp[len(temp)/2:]))
    #else:
    #    peaks.append(2.5)
    
    #Just using this line alone caused massive outliers
    #peaks.append(max(temp[len(temp)/2:]))

#Arrange them in chunks of 20 for the heatmap
grid = []
gridmask = []
for i in range(50):
    grid.append(peaks[50*i:50*(i+1)])
    gridmask.append(maskdata[50*i:50*(i+1)])

numpygrid = np.array(grid)

#BEGIN EXPERIMENTAL CODE
fig,ax = plt.subplots()
#ax.invert_yaxis()
#ax.set_xbound(0,10)
ax.set_xticklabels(xticks)
ax.set_yticklabels(yticks)
#END EXPERIMENTAL CODE


#Legal "cmap" values for color choice are as follows:

# Accent, Accent_r, Blues, Blues_r, BrBG, BrBG_r, BuGn, BuGn_r, BuPu, BuPu_r, CMRmap, CMRmap_r, Dark2, Dark2_r, GnBu, GnBu_r, Greens, Greens_r, Greys, Greys_r, OrRd, OrRd_r, Oranges, Oranges_r, PRGn, PRGn_r, Paired, Paired_r, Pastel1, Pastel1_r, Pastel2, Pastel2_r, PiYG, PiYG_r, PuBu, PuBuGn, PuBuGn_r, PuBu_r, PuOr, PuOr_r, PuRd, PuRd_r, Purples, Purples_r, RdBu, RdBu_r, RdGy, RdGy_r, RdPu, RdPu_r, RdYlBu, RdYlBu_r, RdYlGn, RdYlGn_r, Reds, Reds_r, Set1, Set1_r, Set2, Set2_r, Set3, Set3_r, Spectral, Spectral_r, Wistia, Wistia_r, YlGn, YlGnBu, YlGnBu_r, YlGn_r, YlOrBr, YlOrBr_r, YlOrRd, YlOrRd_r, afmhot, afmhot_r, autumn, autumn_r, binary, binary_r, bone, bone_r, brg, brg_r, bwr, bwr_r, cool, cool_r, coolwarm, coolwarm_r, copper, copper_r, cubehelix, cubehelix_r, flag, flag_r, gist_earth, gist_earth_r, gist_gray, gist_gray_r, gist_heat, gist_heat_r, gist_ncar, gist_ncar_r, gist_rainbow, gist_rainbow_r, gist_stern, gist_stern_r, gist_yarg, gist_yarg_r, gnuplot, gnuplot2, gnuplot2_r, gnuplot_r, gray, gray_r, hot, hot_r, hsv, hsv_r, inferno, inferno_r, jet, jet_r, magma, magma_r, nipy_spectral, nipy_spectral_r, ocean, ocean_r, pink, pink_r, plasma, plasma_r, prism, prism_r, rainbow, rainbow_r, seismic, seismic_r, spectral, spectral_r, spring, spring_r, summer, summer_r, terrain, terrain_r, viridis, viridis_r, winter, winter_r
colormap = matplotlib.colormaps.get_cmap("coolwarm")
colormap.set_bad(color="black")

heatmap = plt.imshow(numpygrid,cmap=colormap,interpolation='nearest')
#masked = np.ma.masked_where(gridmask == 0, gridmask)
#ax.imshow(masked, alpha = 1, cmap = "gray")
#Candidates: coolwarm, rainbow
ax = plt.gca()
cbar = ax.figure.colorbar(heatmap, ax=ax)
cbar.ax.set_ylabel("Peak BAX level", rotation=-90, va="bottom")

plt.xlabel(xlabel)
plt.ylabel(ylabel)

plt.show()
