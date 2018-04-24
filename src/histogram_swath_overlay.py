# -*- coding: utf-8 -*-
"""Test program using Python wrapper for timsdata.dll"""

import sys
import pyopenms
import numpy as np, matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


spacing = 0.001128554344177246 * 1.001
spacing = 0.001128554344177246 * 1.0
spacing = 0.001128554344177246 * 0.999
myrange = 1.8
nbin = myrange / spacing
bins = range(int(nbin+1)) 
xaxis = [i*spacing for i in bins]

output_pdf = sys.argv[1]

def get_hist(e):
    hist = [0 for i in bins]
    #
    fda_values = []
    for s in e:
        ims = list(s.getFloatDataArrays()[0])
        mzs, intens = s.get_peaks()
        for im, mz, inten in zip(ims, mzs, intens):
            hist[ int(im/spacing) ] += inten
    return hist

allhist = []
legends = []
for i in range(0,  9):
    print i
    e = pyopenms.MSExperiment()
    # filename = "/tmp/hela_spl/openswath_tmpfile_%s.mzML" % i
    filename = "/tmp/hela_all_split/openswath_tmpfile_%s.mzML" % i
    pyopenms.MzMLFile().load(filename, e)
    h = get_hist(e)
    allhist.append(h)
    s = e[0]
    plt_title = "SWATH: %s - %s" % (s.getPrecursors()[0].getMZ()-12.5, s.getPrecursors()[0].getMZ()+12.5)
    legends.append(plt_title)


def doPlot(xaxis, hist, plt, label):
    # Plot without empty bins
    xpl = []
    hpl = []
    for x, h in zip(xaxis, hist):
        if h > 500:
            xpl.append(x)
            hpl.append(h)
    plt.plot(xpl, hpl, label=label)

plot1 = plt.figure(1, figsize=(10,5))
plt.clf()
plt.title("All Windows")
plt.xlabel("Drift time (1/K0)")
plt.ylabel("Intensity")
plt.xlim(0.6, 1.7)
for h,l in zip(allhist, legends):
    doPlot(xaxis, h, plt, l)

plt.legend(loc=5, bbox_to_anchor=[1, 0.5],
                   ncol=1, shadow=True, title="Legend", fancybox=True)
# plt.show()

pp = PdfPages(output_pdf)
pp.savefig(plot1)
pp.close()


