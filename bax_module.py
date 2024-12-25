# Copyright (c) 2015 Derrick Sund
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.


# This should not be run without an argument, and that argument should be a
# directory containing output from goldbeter_full_ink4_p53.py as well as
# the output of analysis scripts figuring out what the cell cycle is doing.


from scipy.integrate import odeint
from numpy import arange
import time
import datetime
import os
import sys

dirname = ""
ink4_data = []
p53_data = []
me_data = []
ma_data = []
mb_data = []
p27_data = []
cdc20a_data = []

steady_state = False

if len(sys.argv) == 0:
    sys.stderr.write("Use: python bax_module.py [directory]")
    sys.stderr.write("Generates timecourse data for the Bax apoptosis module,")
    sys.stderr.write("given existing cell cycle data in [directory].")

if len(sys.argv) > 1:
    dirname = sys.argv[1]
    if not os.path.isdir(dirname):
        sys.stderr.write(dirname + ": No such directory found.\n")
        sys.exit(1)
    os.chdir(dirname)
    if not os.path.isfile("Ink4.txt"):
        sys.stderr.write("Ink4.txt not found.  Is this directory correct?")
    for line in open("Ink4.txt",'r').readlines():
        ink4_data.append(float(line.strip().split()[1]))

    if not os.path.isfile("p53.txt"):
        sys.stderr.write("p53.txt not found.  Is this directory correct?")
    for line in open("p53.txt",'r').readlines():
        p53_data.append(float(line.strip().split()[1]))

    if not os.path.isfile("Me.txt"):
        sys.stderr.write("Me.txt not found.  Is this directory correct?")
    for line in open("Me.txt",'r').readlines():
        me_data.append(float(line.strip().split()[1]))

    if not os.path.isfile("Ma.txt"):
        sys.stderr.write("Ma.txt not found.  Is this directory correct?")
    for line in open("Ma.txt",'r').readlines():
        ma_data.append(float(line.strip().split()[1]))

    if not os.path.isfile("Mb.txt"):
        sys.stderr.write("Mb.txt not found.  Is this directory correct?")
    for line in open("Mb.txt",'r').readlines():
        mb_data.append(float(line.strip().split()[1]))

    if not os.path.isfile("p27.txt"):
        sys.stderr.write("p27.txt not found.  Is this directory correct?")
    for line in open("p27.txt",'r').readlines():
        p27_data.append(float(line.strip().split()[1]))


    if not os.path.isfile("Cdc20a.txt"):
        sys.stderr.write("Cdc20a.txt not found.  Is this directory correct?")
    for line in open("Cdc20a.txt",'r').readlines():
        cdc20a_data.append(float(line.strip().split()[1]))
    #What we need from previous steps:
    #Ink4
    #p53
    #cell cycle metadata (ma, p27)

#Check for actively dividing cells
me_slice = me_data[(len(me_data)*4)//5:-1]
ma_slice = ma_data[(len(ma_data)*4)//5:-1]
mb_slice = mb_data[(len(mb_data)*4)//5:-1]

if max(me_slice) < 1.2*min(me_slice) or max(ma_slice) < 1.2*min(ma_slice) or \
   max(mb_slice) < 1.2*min(mb_slice) or max(me_slice) < 0.01 or \
   max(ma_slice) < 0.01 or max(mb_slice) < 0.01:
    steady_state = True
    ss = open("STEADY.txt",'w')
    ss.write("TRUE")
    ss.close


cycling_coefficient = 0
if not steady_state:
    #Cycling coefficient: how distorted the cell cycle is.  Ranges from 0.3 to 1.
    cycling_coefficient = float(open("maxratio.txt").read().strip())


#Constants.  Do not add constants directly to the derivative function; violators
#will have rabid weasels set upon them.
eps = 17

beta_p53bcl2 = 0.2
delta_p53bcl2 = 0.2
alpha_bh3 = 0.2
omega_bh3 = 0.15
alpha_bh3p53 = 0.2
beta_bh3bcl2 = 0.2
delta_bh3bcl2 = 0.2
alpha_bcl2 = 0.2
omega_bcl2 = 0.2
alpha_bcl2p16 = 0.2
beta_baxbcl2 = 0.2
delta_baxbcl2 = 0.2
alpha_bax = 0.2
alpha_baxp53 = 0.2
alpha_baxbh3 = 0.2
omega_bax = 0.2
beta_baxp16 = 0.2
delta_baxp16 = 0.2
p16tot = 1
omega_baxp16 = 0.2
omega_baxbcl2 = 0.2
omega_p53bcl2 = 0.2
omega_bh3bcl2 = 0.2

alpha_proofreading = 0.2
omega_proofreading = 10.1

alpha_stress = 0.1
omega_stress_proofreading = 1
omega_stress = 0.05
alpha_baxstress = 4

y0 = [0.1,0.1,0.1,0.1,0.1,0.1,0.1]


names = []

names.append("BH3")
names.append("BCL2")
names.append("BAX")
names.append("p16-BAX")
names.append("BCL2-BAX")
names.append("p53-BCL2")
names.append("BH3-BCL2")


#y[0] = BH3
#y[1] = BCL2
#y[2] = BAX
#y[3] = p16-BAX
#y[4] = BCL2-BAX
#y[5] = p53-BCL2
#y[6] = BH3-BCL2


#Timecourse data is written every 0.1h, so finding 100t rounded to the
#nearest integer gives us the index we want
def p53(t):
    global p53_data
    index = round(10*t)
    if index >= len(p53_data):
        index = -1
    return p53_data[index]

def ink4(t):
    global ink4_data
    index = round(10*t)
    if index >= len(ink4_data):
        index = -1
    return ink4_data[index]

#The derivative function for the differential equation system.
def func(y,t):
    return [
             (alpha_bh3 - omega_bh3 * y[0] - beta_bh3bcl2 * y[0] * y[1] + delta_bh3bcl2 * y[6]) * eps,
             (alpha_bcl2 - omega_bcl2 * y[1] + alpha_bcl2p16 * (ink4(t) - y[3]) - beta_p53bcl2 * (p53(t) - y[5]) * y[1] + delta_p53bcl2 * y[5] - beta_bh3bcl2 * y[0] * y[1] + delta_bh3bcl2 * y[6] - beta_baxbcl2 * y[1] * y[2] + delta_baxbcl2 * y[4]) * eps,
             (alpha_bax + alpha_baxstress * cycling_coefficient + alpha_baxp53 * (p53(t) - y[5]) + alpha_baxbh3 * y[0] - omega_bax * y[2] - beta_baxbcl2 * y[1] * y[2] + delta_baxbcl2 * y[4] - beta_baxp16 * y[2] * (ink4(t) - y[3]) + delta_baxp16 * y[3]) * eps,
             (beta_baxp16 * y[2] * (ink4(t) - y[3]) - delta_baxp16 * y[3] - omega_baxp16 * y[3]) * eps,
             (beta_baxbcl2 * y[1] * y[2] - delta_baxbcl2 * y[4] - omega_baxbcl2 * y[4]) * eps,
             (beta_p53bcl2 * (p53(t) - y[5]) * y[1] - delta_p53bcl2 * y[5] - omega_p53bcl2 * y[5]) * eps,
             (beta_bh3bcl2 * y[0] * y[1] - delta_bh3bcl2 * y[6] - omega_bh3bcl2 * y[6]) * eps,
           ]
           
           
t = arange(0, 500.0, 0.1)

y = odeint(func, y0, t, ixpr=True)


for i in range(len(y0)):
    writer = open(names[i]+".txt", 'w')
    for j in range(len(t)):
        writer.write(str(t[j]) + " " + str(y[j][i]) + "\n")
    writer.close()
