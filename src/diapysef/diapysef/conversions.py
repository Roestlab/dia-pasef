import sys
import matplotlib
matplotlib.use('Agg')
import pandas as pd
import statsmodels.formula.api as sm
import numpy as np
from statsmodels.graphics.regressionplots import abline_plot
from matplotlib.backends.backend_pdf import PdfPages
from scipy.interpolate import interp1d

def calibrate(data, plot = True, pdf = "rtcalibration"):
    """Fits a linear model to the irts"""
    filename = data.iloc[0,0]
    mod = sm.ols(formula = 'irt ~ rt', data = data)
    res = mod.fit()
    # scatter-plot data
    ax = data.plot(x='rt', y='irt', kind='scatter')
    a = abline_plot(model_results=res, ax=ax)
    a = abline_plot(intercept=0, slope=1, ax=ax)
    # print(res.rsquared)
    a.suptitle("%s \n R2: %s \n R2adj: %s" % (filename, res.rsquared, res.rsquared_adj), fontsize = 10)
    # text(0.9, 0.1,("R2:"), ha='center', va='center', transform=ax.transAxes)
    a.savefig(pdf, format = 'pdf')
    matplotlib.pyplot.close()
    return [filename, res.params.Intercept, res.params.rt]
    # row = np.array([['raw', 'intercept', 'slope'], [data.iloc[0,0], res.params.Intercept, res.params.rt]])
    # return pd.DataFrame(data = row[1:,], columns = row[0,0:])

def outliers(data, tolerance = 7):
    '''
    (dataframe, int) -> list of dataframe & list of outliers indices

    Identify outlier indices for lowess alignment.
    Returns a pandas dataframe wwith outliers removed
    '''
    filename = data.iloc[0,0]
    mod = sm.ols(formula = 'irt ~rt', data = data)
    res = mod.fit()
    absd = abs(res.resid) # get the absolute s.d. for each point from linear regression
    outliers = data[absd > 7].index #get the index of the outliers
    data1 = data.drop(outliers)
    return [data1,outliers]

def get_product_charge(msl, msms):
    return(msl)

def pasef_to_tsv(evidence, msms,
                 irt_file = None,
                 pdfout = "rtcalibration.pdf",
                 im_column = 'Ion mobility index',
                 rt_alignment = 'nonlinear',
                 im_alignment = 'linear'):
    """Converts a mq output to a library taking a best replicate approach."""

    ev = evidence.loc[:, ["id", "Calibrated retention time", "Ion mobility index", im_column]]
    ev = ev.rename(columns = {'id':'Evidence ID'})
    ev = ev.rename(columns = {im_column:'imcol'})
    ms = pd.merge(msms, ev, on = 'Evidence ID')
    ms = ms.dropna(subset=["Calibrated retention time"]) # Some precursors in MQ are not annotated with a RT

    # Replace modifications in the MQ output so that they are OpenMS readable
    ms['Modified sequence'] = ms['Modified sequence'].str.replace('_', '')
    ms['Modified sequence'] = ms['Modified sequence'].str.replace("(ox)","Oxidation")
    ms['Modified sequence'] = ms['Modified sequence'].str.replace("(ph)","Phospho")
    ms['Modified sequence'] = ms['Modified sequence'].str.replace("C","C(Carbamidomethyl)")
    ms['Modified sequence'] = ms['Modified sequence'].str.replace("(ac)","Acetylation")
    ms = ms.rename(columns = {'Modified sequence':'ModifiedPeptideSequence'})

    if irt_file is not None:
        if isinstance(irt_file, str):
            irt = pd.read_table(irt_file)
        elif isinstance(irt_file, pd.DataFrame):
            irt = irt_file
        else:
            print("irt_file must be a path to an irt table or a pd.DataFrame object.")
            sys.exit()

            # Prepare iRT table
        if irt.shape[1] > 2:
            irt = irt.loc[:, ["ModifiedPeptideSequence", "NormalizedRetentionTime", "PrecursorIonMobility","PrecursorCharge"]]
            irt = irt.drop_duplicates()

        irt.columns = ["sequence","irt", "iim", "charge"]
        irt['sequence'] = irt['sequence'].str.replace("\\(UniMod:4\\)","C(Carbamidomethyl)")
        irt['sequence'] = irt['sequence'].str.replace("\\(UniMod:1\\)","Acetylation")
        irt['sequence'] = irt['sequence'].str.replace("\\(UniMod:35\\)","Oxidation")
        irt['sequence'] = irt['sequence'].str.replace("\\(UniMod:21\\)","Phospho")

        # Use the irt peptides present in the MQ output to calibrate
        msms_irt = ms[ms["ModifiedPeptideSequence"].isin(irt["sequence"])]
        ## We chose the best id of each peptide (only one charge state for alignment)
        msms_irt = msms_irt.loc[msms_irt.groupby(["Raw file","ModifiedPeptideSequence"])['PEP'].idxmin()]
        msms_irt = msms_irt.loc[:, ["Raw file","ModifiedPeptideSequence","Calibrated retention time", "imcol", "Charge"]]
        msms_irt.columns = ["raw","sequence","rt","im", "charge"]
        msms_irt = pd.merge(msms_irt, irt, on = ['sequence', 'charge'])
        raw_files = msms_irt.raw.unique()
       
        if rt_alignment is not None:
            if rt_alignment is 'linear':
                print("Aligning retention time linearly...")
                # Generate the iRT calibrators
                calibrators = []
                pp = PdfPages(pdfout)
                for file in raw_files:
                    cal = calibrate(msms_irt[msms_irt.raw == file], pdf = pp)
                    calibrators.append(cal)
                pp.close()

                # Apply rt calibrators
                calibrators = pd.DataFrame(columns = ['Raw file', 'intercept', 'slope'], data = calibrators)
                ms = pd.merge(ms, calibrators, on = 'Raw file')
                ms['irt'] = ms.intercept + ms.slope * ms['Calibrated retention time']
                ms = ms.drop(columns=['intercept', 'slope'])

            elif rt_alignment is 'nonlinear':
                print("Aligning retention time by lowess...")
                import statsmodels.api as smnonlinear
                lowess = smnonlinear.nonparametric.lowess
                # make a df to store the fitted values for merging later
                for file in raw_files:
                    msms_irt_sub = msms_irt[msms_irt.raw == file]
                    msms_irt_sub = msms_irt_sub.loc[:,["rt","irt"]]
                    # does not need to remove outliers
                    # span for the lowess interpolation is 1% of the range
                    delta = (max(msms_irt_sub['rt'] - min(msms_irt_sub['rt'])) * 0.01 )
                    if len(msms_irt_sub) < 100:
                        frac = 1.0 # take the fraction of data to generate the weighted fitting
                    else:
                        frac = 0.1
                    # lowess fits by y/x
                    r = lowess(msms_irt_sub['irt'], msms_irt_sub['rt'], delta=delta, frac=frac)
                    lowess_x = list(zip(*r))[0]
                    lowess_y = list(zip(*r))[1]
                    # create an interpolation function
                    f = interp1d(lowess_x, lowess_y, bounds_error=False)
                    nRT = np.asarray(f(ms.loc[ms['Raw file'] == file ,'Calibrated retention time'].values))

                    # fit linear extrapolation for values outside of approximation
                    # print(sum(np.isnan(nRT)))
                    min_bound = min(lowess_x)
                    max_bound = max(lowess_x)
                    idx = np.asarray(np.where((ms.loc[ms['Raw file'] == file, 'Calibrated retention time'].values < min_bound) | (ms.loc[ms['Raw file'] == file, 'Calibrated retention time'].values > max_bound)))

                    # cannot use pandas.index on ms, it returns the index of the original ms data frame and not th subsetted dataframe

                    lnmod = calibrate(msms_irt_sub)
                    intercept = lnmod[1] # intercept
                    slope = lnmod[2] # slope
                    nRT[idx] = slope * (ms.loc[ms['Raw file'] == file, 'Calibrated retention time'].values[idx]) + intercept
     
                    nrt = []
                    for t in nRT: nrt.append(t)
                    ms.loc[ms['Raw file'] == file, 'irt'] = list(map(str, nrt))
            else:
                print("Only rt_alignment:linear and lowess calibrations are currently implemented")
                ms['irt'] = ms['Calibrated retention time']
        else:
            ms['irt'] = ms['Calibrated retention time']

        if im_alignment is not None:
            if im_alignment is "linear":
                print("Aligning ion mobility")
                # Generate the iIM calibrators
                calibrators = []
                pp = PdfPages("imcalibration.pdf")
                for file in raw_files:
                    msms_irt_sub = msms_irt[msms_irt.raw == file]
                    msms_irt_sub = msms_irt_sub.loc[:,["raw","sequence","im","iim"]]
                    msms_irt_sub.columns = ["raw","sequence","rt","irt"]
                    cal = calibrate(msms_irt_sub, pdf = pp)
                    calibrators.append(cal)
                pp.close()
                # Apply im calibrators
                calibrators = pd.DataFrame(columns = ['Raw file', 'intercept', 'slope'], data = calibrators)
                ms = pd.merge(ms, calibrators, on = 'Raw file')
                ms['iim'] = ms.intercept + ms.slope * ms['imcol']
                ms = ms.drop(columns=['intercept', 'slope'])
            else:
                print("Only im_alignment='linear' is currently implemented")
                ms['iim'] = ms['imcol']

        else:
            ms['iim'] = ms['imcol']
    else:
        ms['irt'] = ms['Calibrated retention time']
        ms['iim'] = ms['imcol']

    # Filter for best identification (lowest PEP)
    ms = ms.loc[ms.groupby(["ModifiedPeptideSequence", "Charge"])['PEP'].idxmin()]
    # Filter decoy peptides
    ms = ms[ms['Reverse'].isna()]

    # Shape data for AssayGenerator
    msl = ms.loc[:, ["id","m/z","Masses","Charge","irt", "iim", "imcol", "Intensities","Sequence","ModifiedPeptideSequence","Proteins"]]
    masses = msl['Masses'].str.split(';', expand = True).stack().str.strip().reset_index(level = 1, drop=True)
    intensities = msl['Intensities'].str.split(';', expand = True).stack().str.strip().reset_index(level = 1, drop=True)
    df1 = pd.concat([masses, intensities], axis = 1, keys = ['Masses', 'Intensities'])
    msl2 = msl.drop(['Masses', 'Intensities'], axis = 1)
    msl2 = msl2.join(df1).reset_index(drop = True)
    msl2.columns = ["transition_group_id","PrecursorMz","PrecursorCharge","iRT", "Im_recalibrated", "PrecursorIonMobility", "PeptideSequence","FullUniModPeptideName","ProteinName", "ProductMz", "LibraryIntensity"]
    msl2 = get_product_charge(msl2, msms)

    # Reorder the columns as they are in libraries generated with the OSW assay generator
    # msl2 = msl2[["transition_group_id","PrecursorMz","ProductMz","PrecursorCharge","Tr_recalibrated","LibraryIntensity","PeptideSequence","FullUniModPeptideName","ProteinName", "PrecursorIonMobility"]]
    msl2['decoy'] = 0
    msl2['transition_name'] = ['_'.join(str(i) for i in z) for z in zip(msl2.transition_group_id,msl2.PrecursorMz,msl2.ProductMz)]

    msl2=msl2.drop_duplicates()
    # Write output
    # msl2.to_csv("pasefLib.tsv", sep = "\t", index=False)

    return(msl2)

# @TODO call the OpenSwathAssayGenerator with pyopenms
# def run_assay_generator(paseflib, window_scheme):
#     import pyopenms

