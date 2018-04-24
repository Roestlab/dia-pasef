
import sys
import numpy as np
import pandas as pd
import seaborn as sns
import timsdata
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from scipy.stats import kde

analysis_dir = '20180320_AnBr_SA_diaPASEF_200ng_HeLa_Rost_Method_4_a_01_A1_01_2143.d'
if sys.version_info.major == 2:
    analysis_dir = unicode(analysis_dir)

td = timsdata.TimsData(analysis_dir)
conn = td.conn

# Get total frame count:
q = conn.execute("SELECT COUNT(*) FROM Frames")
row = q.fetchone()
N = row[0]
print("Analysis has {0} frames.".format(N))

def get_cycle_length(conn, pasef_ms2_id = 8):
    firstcycle = True
    initial_ms1 = True
    cycle_length = 0
    row = 1
    while (firstcycle):
        mstype = conn.execute("SELECT MsMsType FROM Frames LIMIT {0},1".format(row)).fetchone()
        row += 1
        if mstype[0] == pasef_ms2_id:
            cycle_length += 1
            initial_ms1 = False
        elif(initial_ms1 != True):
            firstcycle = False
    return cycle_length

def get_windows_per_frame(conn, pasef_ms2_id = 8):
    q = conn.execute("SELECT * FROM Frames INNER JOIN PasefFrameMsMsInfo ON Frames.Id=PasefFrameMsMsInfo.Frame WHERE MsMsType=%d AND Frame = (SELECT MIN(Frame) FROM PasefFrameMsMsInfo as fr)" % pasef_ms2_id)
    windows = q.fetchall()
    return len(windows)

def get_conversion_func(conn):
    q = conn.execute("SELECT * FROM TimsCalibration")
    calib = q.fetchone()
    def convert_im(im):
        im = np.array(im)
        return(1/(calib[8]+calib[9]/(calib[4]+((calib[5]-calib[4])/calib[3])*(im-calib[6]-calib[2]))))
    # x = calib
    # Mobility[1/k0] = 1/(c6+c7/(c2+((c3-c2)/c1)*(scanno-c4-c0)))
    return convert_im

convert_scan_num = get_conversion_func(conn)
convert_scan_num([ 100, 200 ])
wpf = get_windows_per_frame(conn, 8)
cl = get_cycle_length(conn)

def get_windows(timsData):
    """Extracts the window scheme from the first cycle of a tims file"""
    conn = timsData.conn
    pasef_ms2_id = 8 # diaPASEF ms2 scans are denoted by 8 instead of 2
    cycle_length = get_cycle_length(conn, pasef_ms2_id)
    wpf = get_windows_per_frame(conn, pasef_ms2_id)
    q = timsData.conn.execute("SELECT * FROM Frames INNER JOIN PasefFrameMsMsInfo ON Frames.Id=PasefFrameMsMsInfo.Frame WHERE MsMsType=%d LIMIT %d" % (pasef_ms2_id, cycle_length*wpf))
    frames = q.fetchall()
    colnames = [description[0] for description in q.description]
    resframe = pd.DataFrame(data = frames, columns = colnames)
    return resframe

w = get_windows(td)
a = np.array(w['ScanNumBegin'])
w['ImBegin'] = convert_scan_num(np.array(w['ScanNumBegin']))
w['ImEnd'] = convert_scan_num(np.array(w['ScanNumEnd']))

# sns_plot1 = sns.jointplot(x=w["IsolationMz"], y=w["ImBegin"], kind='scatter')
# a = w['ImEnd'][2]
iris = sns.load_dataset("iris")
f, ax = plt.subplots(figsize = (5,6))
i =1
for p in [
    patches.Rectangle(
        (w['IsolationMz'][i]-w['IsolationWidth'][i]/2, w['ImBegin'][i]),
        w['IsolationWidth'][i],
        w['ImEnd'][i] - w['ImBegin'][i]) for i in range(len(w))
]:
    ax.add_patch(p)
# ax.add_patch(patches.Rectangle((0.400,0.500), 0.5,0.8))
ax.set(xlim = (400, 1200),
       ylim = (0.5,2))
# Add rectangle
# sns_plot1.add_patch(
# patches.Rectangle(
# (20, 25), # (x,y)
# 50, # width
# 50, # height
# # You can add rotation as well with 'angle'
# alpha=0.3, facecolor="red", edgecolor="black", linewidth=3, linestyle='solid'
# )
# )
# sns_plot1.savefig("output.png")
plt.show()

msone = pd.read_table("allPeptides.txt")
msone = msone[msone.Charge == 2]
msone = msone[['m/z', 'Ion mobility index', 'Intensity']]
msone.columns = ["mz", "scannr", "intensity"]
nr_bins = 2000
msone["binmz"] = pd.cut(msone.mz, nr_bins)
msone["im"] = convert_scan_num(msone.scannr)
msone["binim"] = pd.cut(msone.im, nr_bins)
# f, ax = plt.subplots(figsize = (8,6))
# ax.scatter(msone.mz, msone.im)
# ax.set(xlim = (400, 1200),
#        ylim = (0.5,2))
# plt.show()

# sns.jointplot(x=msone['mz'], y = msone['im'], kind = 'kde')

# Evaluate a gaussian kde on a regular grid of nbins x nbins over data extents
# nbins=200
# k = kde.gaussian_kde([msone.mz, msone.im])
# xi, yi = np.mgrid[msone.mz.min():msone.mz.max():nbins*1j, msone.im.min():msone.im.max():nbins*1j]
# zi = k(np.vstack([xi.flatten(), yi.flatten()]))
 
# # Make the plot
# plt.pcolormesh(xi, yi, zi.reshape(xi.shape))
plt.hist2d(msone.mz, msone.im)
plt.show()


