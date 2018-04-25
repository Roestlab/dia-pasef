import pandas as pd
import os.path

def calculate_coverage (windows, precursor_map = None, count_sc = False):
    """Calculates how many precursors are covered by the window setup."""
    if precursor_map is None:
        mq = pd.read_pickle(os.path.join(os.path.dirname(__file__), 'data/all_peptides.pk'))
    else:
        mq = precursor_map
    if not count_sc:
        mq = mq[mq.Charge > 1]
    w = windows
    w['mzstart'] = w['IsolationMz'] - w['IsolationWidth']/2
    w['mzend'] = w['IsolationMz'] + w['IsolationWidth']/2
    mq['is_covered'] = False
    for i in range(len(windows)):
        mq.loc[(mq['m/z'] < w.loc[i,'mzend']) & (mq['m/z'] > w.loc[i, 'mzstart']) & (mq['IonMobilityIndexK0'] > w.loc[i, 'IMend']) & (mq['IonMobilityIndexK0'] < w.loc[i, 'IMstart']),'is_covered'] = True
    # cov = sum(mq.is_covered)/len(mq)
    return mq
