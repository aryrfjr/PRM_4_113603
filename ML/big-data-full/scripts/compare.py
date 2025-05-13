#!/usr/bin/python

import sys
import numpy as np
import pickle as pkl
import matplotlib.pyplot as plt

import math

ftl = open("/home/aryjr/UFSCar/MG-NMR/ML/big-data-full/Zr45Cu45Al10-PBSSDB/Al-Al-SOAPS.vec","rb")
sss = pkl.load(ftl)
ftl.close()
ftl = open("/home/aryjr/UFSCar/MG-NMR/ML/big-data-full/Zr45Cu45Al10-PBMSDB/Al-Al-SOAPS.vec","rb")
mss = pkl.load(ftl)
ftl.close()
ss = sss["1-0-99"]
ms = mss["1-0-99"]
print np.dot(ss,ms)

