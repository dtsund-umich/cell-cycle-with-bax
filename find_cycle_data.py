import os
import sys
import subprocess
import bisect

if len(sys.argv) < 2:
    print("Try that again but with an argument next time.")
    sys.exit(1)

#Autocomplete might give us a slash in the directory name; this leads
#to confusing behavior if we don't strip it out.
model = sys.argv[1].split("/")[0]

#Step 2: Load Cdc20a, p27, and Ma to determine phase transition points

#G1/S-G2 transition: peak p27?
#S-G2/G2-M transition: peak Ma?
#G2-M/G1 transition: peak Cdc20a?

p27lines = open(model+"/p27.txt").readlines()
malines = open(model+"/Ma.txt").readlines()
cdc20alines = open(model+"/Cdc20a.txt").readlines()

p27data = []
madata = []
cdc20adata = []

for i in range(len(p27lines)):
    p27data.append(float(p27lines[i].split()[1]))
    madata.append(float(malines[i].split()[1]))
    cdc20adata.append(float(cdc20alines[i].split()[1]))
    

#For now, I will assume one peak per cycle.
#This will break in very visible ways if the assumption doesn't hold, and I'll go back and fix it then if it ever comes up.
g1_sg2_trans = []
sg2_g2m_trans = []
g2m_g1_trans = []


#Start in g1.  Doesn't really matter, initial transient behavior will get smoothed out.
g1_data = [1]
sg2_data = [0]
g2m_data = [0]

#Note what the peak of each of the three timecourses is after initial transients have been ironed out.
#Peaks only actually occur when we're near this level, some spurious "peaks" can occur elsewhere
p27peak = max(p27data[len(p27data)//2:])
mapeak = max(madata[len(madata)//2:])
cdc20apeak = max(cdc20adata[len(cdc20adata)//2:])

if p27peak < 0.01 or mapeak < 0.01 or cdc20apeak < 0.01:
    open(model+"/period.txt", 'w').write("0")
    open(model+"/maxratio.txt", 'w').write("0")
    sys.exit()

current_phase = "g1"

#Minor assumption: all three loaded files are the same length.  If this doesn't hold, something is very seriously wrong...
for i in range(1,len(p27lines)-1):
    #if p27data[i] > p27data[i-1] and p27data[i] > p27data[i+1] and p27data[i]/p27peak > 0.9:
    if p27data[i] == max(p27data[max([0,i-10]):min([i+10,len(p27lines)-1])]) and p27data[i]/p27peak > 0.9 and current_phase == "g1":
        #Transition to S/G2
        current_phase = "sg2"
        g1_sg2_trans.append(float(p27lines[i].split()[0]))
    #if madata[i] > madata[i-1] and madata[i] > madata[i+1] and madata[i]/mapeak > 0.9:
    if madata[i] == max(madata[max([0,i-10]):min([i+10,len(malines)-1])]) and madata[i]/mapeak > 0.9 and current_phase == "sg2":
        #Transition to G2/M
        current_phase = "g2m"
        sg2_g2m_trans.append(float(malines[i].split()[0]))
    #if cdc20adata[i] > cdc20adata[i-1] and cdc20adata[i] > cdc20adata[i+1] and cdc20adata[i]/cdc20apeak > 0.9:
    if cdc20adata[i] == max(cdc20adata[max([0,i-10]):min([i+10,len(cdc20alines)-1])]) and cdc20adata[i]/cdc20apeak > 0.9 and current_phase == "g2m":
        #Transition to G1
        current_phase = "g1"
        g2m_g1_trans.append(float(cdc20alines[i].split()[0]))
    
    
    if current_phase == "g1":
        g1_data.append(1)
        sg2_data.append(0)
        g2m_data.append(0)
    if current_phase == "sg2":
        g1_data.append(0)
        sg2_data.append(1)
        g2m_data.append(0)
    if current_phase == "g2m":
        g1_data.append(0)
        sg2_data.append(0)
        g2m_data.append(1)
        

if current_phase == "g1":
    g1_data.append(1)
    sg2_data.append(0)
    g2m_data.append(0)
if current_phase == "sg2":
    g1_data.append(0)
    sg2_data.append(1)
    g2m_data.append(0)
if current_phase == "g2m":
    g1_data.append(0)
    sg2_data.append(0)
    g2m_data.append(1)


period = 0
g1_length = 0
sg2_length = 0
g2m_length = 0
maxratio = 0

if len(g1_sg2_trans) > 3:
    #if this isn't true, we're probably in a steady-state solution or one that should be considered such
    period = g1_sg2_trans[-1] - g1_sg2_trans[-2]
    #We don't know which phase we ended in, so we have to do these checks to be sure we're getting the
    #right transitions
    if sg2_g2m_trans[-2] > g1_sg2_trans[-2]:
        sg2_length = sg2_g2m_trans[-2] - g1_sg2_trans[-2]
    else:
        sg2_length = sg2_g2m_trans[-1] - g1_sg2_trans[-2]
    if g2m_g1_trans[-2] > sg2_g2m_trans[-2]:
        g2m_length = g2m_g1_trans[-2] - sg2_g2m_trans[-2]
    else:
        g2m_length = g2m_g1_trans[-1] > sg2_g2m_trans[-2]
    if g1_sg2_trans[-2] > g2m_g1_trans[-2]:
        g1_length = g1_sg2_trans[-2] - g2m_g1_trans[-2]
    else:
        g1_length = g1_sg2_trans[-1] - g2m_g1_trans[-2]

    maxratio = max([g1_length,sg2_length,g2m_length])/period
    
    if maxratio > 1:
        print("Error in " + model + "; maxratio is " + str(maxratio) + ", greater than 1")
    if maxratio < 0:
        print("Error in " + model + "; maxratio is " + str(maxratio) + ", greater than 1")



#Write all the potentially relevant data.
g1_writer = open(model+"/g1.txt", 'w')
sg2_writer = open(model+"/sg2.txt", 'w')
g2m_writer = open(model+"/g2m.txt", 'w')
for i in range(len(g1_data)):
    g1_writer.write(p27lines[i].split()[0] + " " + str(g1_data[i]) + "\n")
    sg2_writer.write(p27lines[i].split()[0] + " " + str(sg2_data[i]) + "\n")
    g2m_writer.write(p27lines[i].split()[0] + " " + str(g2m_data[i]) + "\n")

g1_writer.close()
sg2_writer.close()
g2m_writer.close()

open(model+"/period.txt", 'w').write(str(period))
open(model+"/maxratio.txt", 'w').write(str(maxratio))



#X-range for rectangle: 92 up to 704
#Y-range: let's try 60 to 90 for starters


"""
#Step 3: Use imagemagick to mangle original image into finalized image
transtime = g1_sg2_trans[bisect.bisect(g1_sg2_trans, 400)-1]
curphase = "sg2"
if sg2_g2m_trans[bisect.bisect(sg2_g2m_trans, 400)-1] > transtime:
    curphase = "g2m"
    transtime = sg2_g2m_trans[bisect.bisect(sg2_g2m_trans, 400)-1]
if g2m_g1_trans[bisect.bisect(g2m_g1_trans, 400)-1] > transtime:
    curphase = "g1"
    transtime = g2m_g1_trans[bisect.bisect(g2m_g1_trans, 400)-1]

transtime = 400

while True:
    if curphase == "sg2":
        newindex = bisect.bisect(sg2_g2m_trans,transtime)
        endpoint = 0
        if newindex == len(sg2_g2m_trans):
            endpoint = 500
        else:
            endpoint = sg2_g2m_trans[newindex]
        #draw rectangle
        xstart = str((transtime - 400) * 6.13 + 91)
        xend = str((endpoint - 400) * 6.13 + 91)
        command = "convert " + image + " -fill DarkSlateGray2 -stroke black -draw \'rectangle " + xstart + " 60 " + xend + " 90\' " + image
        os.system(command)
        if endpoint == 500:
            break
        transtime = endpoint
        curphase = "g2m"


    elif curphase == "g2m":
        newindex = bisect.bisect(g2m_g1_trans,transtime)
        endpoint = 0
        if newindex == len(g2m_g1_trans):
            endpoint = 500
        else:
            endpoint = g2m_g1_trans[newindex]
        #draw rectangle
        xstart = str((transtime - 400) * 6.13 + 91)
        xend = str((endpoint - 400) * 6.13 + 91)
        command = "convert " + image + " -fill firebrick3 -stroke black -draw \'rectangle " + xstart + " 60 " + xend + " 90\' " + image
        os.system(command)
        if endpoint == 500:
            break
        transtime = endpoint
        curphase = "g1"


    elif curphase == "g1":
        newindex = bisect.bisect(g1_sg2_trans,transtime)
        endpoint = 0
        if newindex == len(g1_sg2_trans):
            endpoint = 500
        else:
            endpoint = g1_sg2_trans[newindex]
        #draw rectangle
        xstart = str((transtime - 400) * 6.13 + 91)
        xend = str((endpoint - 400) * 6.13 + 91)
        command = "convert " + image + " -fill LightGoldenrod -stroke black -draw \'rectangle " + xstart + " 60 " + xend + " 90\' " + image
        os.system(command)
        if endpoint == 500:
            break
        transtime = endpoint
        curphase = "sg2"

    else:
        print("FLAGRANT ERROR")
        break
"""

