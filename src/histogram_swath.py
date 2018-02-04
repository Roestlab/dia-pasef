# -*- coding: utf-8 -*-
"""Test program using Python wrapper for timsdata.dll"""

import sys
import pyopenms
import numpy as np, matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

e = pyopenms.MSExperiment()
filename = "/tmp/hela_spl/openswath_tmpfile_3.mzML"
filename = sys.argv[1]
output_pdf = sys.argv[2]
pyopenms.MzMLFile().load(filename, e)

spacing = 0.001128554344177246 * 1.001
spacing = 0.001128554344177246 * 1.0
spacing = 0.001128554344177246 * 0.999
# spacing = 0.0025
myrange = 1.8
nbin = myrange / spacing
bins = range(int(nbin+1)) 
hist = [0 for i in bins]
xaxis = [i*spacing for i in bins]
hist2d = [ [0 for k in range(2000)] for i in bins]

fda_values = []
for s in e:
    ims = list(s.getFloatDataArrays()[0])
    mzs, intens = s.get_peaks()
    for im, mz, inten in zip(ims, mzs, intens):
        hist[ int(im/spacing) ] += inten
        hist2d[ int(im/spacing) ][int(mz)] += inten

plt_title = "SWATH: %s - %s" % (s.getPrecursors()[0].getMZ()-12.5, s.getPrecursors()[0].getMZ()+12.5)

q = np.array(hist2d)
q = np.log(q+1)
#plot1, ax = plt.figure(1, figsize=(10,5))
plot1, ax1 = plt.subplots()
#plot1.clf()
# ax1.set_xticks(ticks)
# ax1.set_xticklabels(ticklabels)
# ax.set_xticks([0, 100])
# ax.set_xticklabels(["a", "b"])
ax1.set_yticks( range(0, 1600, 350) )
# ax1.set_yticklabels(['', 0,10,20,30,40])
ax1.set_yticklabels( ["{:10.4f}".format(k*spacing) for k in range(0, 1600, 350)] )
ax1.imshow(q, cmap='hot', interpolation='nearest')
plt.ylabel("Drift time (1/K0)")
plt.xlabel("m/z")
#plt.show()
plt.title(plt_title)
# pp = PdfPages("/tmp/a.pdf")
pp = PdfPages(output_pdf + "_heat.pdf")
pp.savefig(plot1)
pp.close()

# plot1 = plt.figure(1, figsize=(10,5))
# plt.clf()
# plt.plot(xaxis, hist)
# plt.xlim(0.6, 1.7)
# pp = PdfPages("/tmp/b.pdf")
# pp.savefig(plot1)
# pp.close()

xpl = []
hpl = []
for x, h in zip(xaxis, hist):
    if h > 500:
        xpl.append(x)
        hpl.append(h)

# plt.figure(1, figsize=(10,5))
# plt.clf()
# # plt.stem(xaxis, hist)
# plt.plot(xaxis, hist)
# plt.xlim(0.6, 1.7)
# plt.show()


plot1 = plt.figure(1, figsize=(10,5))
plt.clf()
plt.title(plt_title)
plt.xlabel("Drift time (1/K0)")
plt.ylabel("Intensity")
plt.xlim(0.6, 1.7)
plt.plot(xpl, hpl)
#plt.plot(xaxis, hist)

pp = PdfPages(output_pdf + "_hist.pdf")
pp.savefig(plot1)
pp.close()


