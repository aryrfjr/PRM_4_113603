#!/usr/bin/python
#

# libraries
import sys
from numpy import array
import matplotlib.pyplot as plt
import matplotlib as mpl

# https://matplotlib.org/2.0.2/examples/color/colormaps_reference.html
norm = mpl.colors.Normalize(vmin=0.0, vmax=2.0)

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2) # , figsize=(11, 3))

# loading data 1
data = []
fileobj = open("1_0_G1Cu_xy.obd")
lines = fileobj.readlines()
data = array([[float(col) for col in line.split()]
            for line in lines[1:len(lines)]])
ax1.matshow(data, cmap=plt.cm.BuPu, norm=norm, aspect='auto')
ax1.set(title='Cu G1')
plt.setp(ax1.get_xticklabels(), visible=False)
plt.setp(ax1.get_yticklabels(), visible=False)

# loading data 2
data = []
fileobj = open("1_0_G1Zr_xy.obd")
lines = fileobj.readlines()
data = array([[float(col) for col in line.split()]
            for line in lines[1:len(lines)]])
ax2.matshow(data, cmap=plt.cm.BuPu, norm=norm, aspect='auto')
ax2.set(title='Zr G1')
plt.setp(ax2.get_xticklabels(), visible=False)
plt.setp(ax2.get_yticklabels(), visible=False)

# loading data 3
data = []
fileobj = open("1_0_G2Zr_xy.obd")
lines = fileobj.readlines()
data = array([[float(col) for col in line.split()]
            for line in lines[1:len(lines)]])
ax3.matshow(data, cmap=plt.cm.BuPu, norm=norm, aspect='auto')
ax3.set(title='Zr G2')
plt.setp(ax3.get_xticklabels(), visible=False)
plt.setp(ax3.get_yticklabels(), visible=False)

# loading data 4
data = []
fileobj = open("1_0_G2Al_xy.obd")
lines = fileobj.readlines()
data = array([[float(col) for col in line.split()]
            for line in lines[1:len(lines)]])
ax4.matshow(data, cmap=plt.cm.BuPu, norm=norm, aspect='auto')
ax4.set(title='Al G2')
plt.setp(ax4.get_xticklabels(), visible=False)
plt.setp(ax4.get_yticklabels(), visible=False)

#plt.colorbar(im)

plt.show()

