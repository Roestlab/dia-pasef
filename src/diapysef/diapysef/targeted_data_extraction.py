import sys
import os

# MS data modules
import pyopenms as po
from .util import setCompressionOptions, method_timer, code_block_timer, setup_logger, check_sqlite_table, check_im_array, type_cast_value


# Logging and performance modules
import traceback
import logging
from tqdm import tqdm

# Modules for data
import pandas as pd
import numpy as np
import itertools
import sqlite3 as sql3
import pickle as pkl
from joblib import Parallel, delayed, wrap_non_picklable_objects

def generate_coordinates(file, outfile=None, run_id=None, target_peptides=None, m_score=0.05, use_transition_peptide_mapping=False, use_only_detecting_transitions=True, verbose=0):
    """
    Generate a dictionary of target coordinates to extract from Raw diaPASEF mzML data

    Params:
      file: a file to generate coordinates from.

    Returns:
      a pickle file containing a dictionary of target coordinates
    """
    # For internal debugging
    if False:
        file = "/media/justincsing/ExtraDrive1/Documents2/Roest_Lab/Github/PTMs_Project/synthetic_pool_timstoff/results/20220824_single_mzml_im_cal_linear/osw/merged.osw"
        target_peptides = ["T(UniMod:21)ELISVSEVHPSR", "TELIS(UniMod:21)VSEVHPSR", "TELISVS(UniMod:21)EVHPSR", "TELISVSEVHPS(UniMod:21)R", "LGDLNY(UniMod:21)LIYVFPDRPK", "LGDLNYLIY(UniMod:21)VFPDRPK",
                           "Y(UniMod:21)VC(UniMod:4)EGPSHGGLPGASSEK", "YVC(UniMod:4)EGPS(UniMod:21)HGGLPGASSEK", "YVC(UniMod:4)EGPSHGGLPGAS(UniMod:21)SEK", "YVC(UniMod:4)EGPSHGGLPGASS(UniMod:21)EK"]
        run_id = 8174980892860667876

    if file.lower().endswith("osw"):
        if verbose == 10:
            logging.debug(
                f"Connecting to OSW identifications database: {file}")
        con = sql3.connect(file)

        if use_transition_peptide_mapping:
            use_transition_map = "INNER JOIN TRANSITION_PEPTIDE_MAPPING ON TRANSITION_PEPTIDE_MAPPING.PEPTIDE_ID = PEPTIDE.ID"
            join_on_transition_table_id = "TRANSITION_PEPTIDE_MAPPING.TRANSITION_ID"
        else:
            use_transition_map = "INNER JOIN TRANSITION_PRECURSOR_MAPPING ON TRANSITION_PRECURSOR_MAPPING.PRECURSOR_ID = PRECURSOR.ID"
            join_on_transition_table_id = "TRANSITION_PRECURSOR_MAPPING.TRANSITION_ID"

        if use_only_detecting_transitions:
            detecting_col = '1'
        else:
            detecting_col = '1, 0'

        if check_sqlite_table(con, "SCORE_IPF"):
            join_on_score_table = "INNER JOIN (SELECT * FROM SCORE_IPF) AS SCORE_TABLE ON SCORE_TABLE.FEATURE_ID = FEATURE.ID"
        elif check_sqlite_table(con, "SCORE_MS2"):
            # Restrict to top peak group rank to generate coordinates for. We also filter based on MS2 QVALUE to reduce the number of peptides to generate coordinates for
            join_on_score_table = "INNER JOIN (SELECT * FROM SCORE_MS2 WHERE RANK ==1 AND QVALUE < %s) AS SCORE_TABLE ON SCORE_TABLE.FEATURE_ID = FEATURE.ID" % (m_score)

        if target_peptides is not None:

            if os.path.isfile(target_peptides):
                logging.info(
                    f"Reading target peptides from file: {target_peptides}")
                with open(target_peptides, 'r') as file:
                    data = file.read().replace('\n', '')
                target_peptides = [peptide.strip()
                                   for peptide in data.split(',')]
            elif type(target_peptides) == str:
                # If a string list is passed to CLI, then convert to literal python list of strings
                target_peptides = type_cast_value(target_peptides)

            logging.info(
                f"Generating coordinates for the following peptides: {target_peptides}")

            sql_query = '''
      SELECT 
      PEPTIDE.MODIFIED_SEQUENCE AS peptide,
      PRECURSOR.PRECURSOR_MZ AS precursor_mz,
      PRECURSOR.CHARGE AS charge,
      TRANSITION.PRODUCT_MZ AS product_mz,
      TRANSITION.CHARGE AS product_charge,
      TRANSITION.ANNOTATION AS product_annotation,
      TRANSITION.DETECTING AS product_detecting,
      FEATURE.EXP_RT AS rt_apex,
      FEATURE.LEFT_WIDTH AS left_width,
      FEATURE.RIGHT_WIDTH AS right_width,
      FEATURE.EXP_IM AS im_apex,
      SCORE_TABLE.QVALUE as qvalue
      FROM PRECURSOR
      INNER JOIN PRECURSOR_PEPTIDE_MAPPING ON PRECURSOR_PEPTIDE_MAPPING.PRECURSOR_ID = PRECURSOR.ID
      INNER JOIN (
        SELECT * 
        FROM PEPTIDE
        WHERE PEPTIDE.MODIFIED_SEQUENCE IN ("%s")
        ) AS PEPTIDE ON PEPTIDE.ID = PRECURSOR_PEPTIDE_MAPPING.PEPTIDE_ID
      %s 
      INNER JOIN (
        SELECT *
        FROM TRANSITION
        WHERE TRANSITION.DETECTING in (%s)
      ) AS TRANSITION ON TRANSITION.ID = %s
      INNER JOIN (
        SELECT *
        FROM FEATURE
        WHERE FEATURE.RUN_ID=%s
      ) AS FEATURE ON FEATURE.PRECURSOR_ID = PRECURSOR.ID
      %s
      ''' % ('","'.join(target_peptides), use_transition_map, detecting_col, join_on_transition_table_id, run_id, join_on_score_table)

            if verbose == 10:
                logging.debug(f"Injecting SQL Query:\n{sql_query}")

            data = pd.read_sql_query(sql_query, con)

            data['group_id'] = data['peptide'] + \
                '_' + data['charge'].astype(str)

            # Only keep the best scoring feature per group_id
            data = data.loc[data.groupby('group_id')['qvalue'].transform(
                'min').eq(data['qvalue'])].reset_index(drop=True)

            # Group and aggregate product m/z into a list
            agg_product_mz_df = data.groupby(['group_id'])['product_mz'].apply(
                list).to_frame().reset_index()
            data = pd.merge(data.drop('product_mz', axis=1),
                            agg_product_mz_df, on=['group_id'])
            # Group and aggregate product charge into a list
            agg_product_charge_df = data.groupby(['group_id'])['product_charge'].apply(
                list).to_frame().reset_index()
            data = pd.merge(data.drop('product_charge', axis=1),
                            agg_product_charge_df, on=['group_id'])
            # Group and aggregate product annotation into a list
            agg_product_annotation_df = data.groupby(['group_id'])['product_annotation'].apply(
                list).to_frame().reset_index()
            data = pd.merge(data.drop('product_annotation', axis=1),
                            agg_product_annotation_df, on=['group_id'])
            # Group and aggregate product detecting into a list
            agg_product_detecting_df = data.groupby(['group_id'])['product_detecting'].apply(
                list).to_frame().reset_index()
            data = pd.merge(data.drop('product_detecting', axis=1),
                            agg_product_detecting_df, on=['group_id'])

            # Aggegate left and right bounaries into a list
            data['rt_boundaries'] = data[[
                'left_width', 'right_width']].values.tolist()
            data.drop(['left_width', 'right_width'], axis=1, inplace=True)

            # Melt dataframe to long format
            data_long = data.melt(id_vars=['group_id'])

            # Generate nested dictionary of targeted peptide coordinates
            peptide_coordinates_dict = data_long.groupby(
                'group_id')[['variable', 'value']].apply(lambda x: dict(x.to_numpy())).to_dict()

        # Close connection
        con.close()

    if outfile is not None:
        logging.info(
            f"Writing coordinates dictionary to pickle file {outfile}")
        with open(outfile, "wb") as output_file:
            pkl.dump(peptide_coordinates_dict, file=output_file)
    else:
        return peptide_coordinates_dict


class data_io():
    """
    Class for data input and output operations
    """

    def __init__(self, mzml_file, mz_tol=20, rt_window=50, im_window=0.06, mslevel=[1], verbose=0, log_file=None):

        # Initialise logger
        if log_file is not None:
            setup_logger(log_file, verbose)

        self.mzml_file = mzml_file
        self.mz_tol = mz_tol
        self.rt_window = rt_window
        self.im_window = im_window
        self.mslevel = mslevel
        self.verbose = verbose
        self.log_file = log_file
        self.consumer = None

    @method_timer
    def load_data(self, is_filtered=False):
        """
        Method to load data from an mzML file as an on disc experiment for memory efficiency and meta data access without loading full data

        Params:
          is_filtered: (boolean) is the input mzML file a filtered mzML file that is already processed for targeted coordinate extraction
        """

        if is_filtered:
            with code_block_timer('Creating MSExperiment...', logging.debug):
                exp = po.MSExperiment()

            with code_block_timer(f'Loading {self.mzml_file} file...', logging.info):
                po.MzMLFile().load(self.mzml_file, exp)

            logging.info(
                f"There are {exp.getNrSpectra()} spectra and {exp.getNrChromatograms()} chromatograms.")
            logging.info(
                f"There are {sum([spec.getMSLevel()==1  for spec in exp.getSpectra()])} MS1 spectra and {sum([spec.getMSLevel()==2  for spec in exp.getSpectra()])} MS2 spectra.")

            self.filtered = exp
        else:
            # TODO: For some reason the OnDiscExperiment doesn't read StringMetaData, which s currecntly used to store filtered IM data
            with code_block_timer('Creating OnDiscExperiment...', logging.debug):
                exp = po.OnDiscMSExperiment()

            with code_block_timer(f'Opening {self.mzml_file} file...', logging.info):
                exp.openFile(self.mzml_file)

            with code_block_timer(f'Extracting meta data...', logging.debug):
                meta_data = exp.getMetaData()

            logging.info(
                f"There are {meta_data.getNrSpectra()} spectra and {exp.getNrChromatograms()} chromatograms.")
            logging.info(
                f"There are {sum([spec.getMSLevel()==1  for spec in meta_data.getSpectra()])} MS1 spectra and {sum([spec.getMSLevel()==2  for spec in meta_data.getSpectra()])} MS2 spectra.")

            self.exp = exp
            self.meta_data = meta_data

    def get_consumer(self, output_fname):
        """
        Get a consumer for mzML or sqMass data writing
        Adpated from diapysef
        """
        # Store output
        if output_fname.lower().endswith("mzml"):
            consumer = po.PlainMSDataWritingConsumer(output_fname)

            # TODO: I've commented out compression, because it seems to make the quality of data bas when uncompressing
            # # Compress output
            # try:
            #     opt = consumer.getOptions()
            #     setCompressionOptions(opt)
            #     consumer.setOptions(opt)
            # except Exception as e:
            #     print(e)
            #     print(
            #         "Your version of pyOpenMS does not support any compression, your files may get rather large")
            #     pass

        elif output_fname.lower().endswith("sqmass"):
            # TODO: For some reason the Sql consumer doesn't work if output file is of type sqmass
            consumer = po.MSDataSqlConsumer(output_fname)

        else:
            raise Exception("Supported filenames: mzML and sqMass.")

        self.consumer = consumer

    @method_timer
    def get_filtered_df(self, ms_level=[1]):
        """
        Get filtered spectrum data into a pandas dataframe format

        Params:
          self: (object) self object that contains filtered data
          ms_level: (list) list of ms level data to write out to file
        """
        results_df = pd.DataFrame()
        for k in tqdm(range(self.filtered.getNrSpectra())):
            spec = self.filtered.getSpectrum(k)
            if spec.getMSLevel() in ms_level:
                mz, intensity = spec.get_peaks()
                rt = np.full([mz.shape[0]], spec.getRT(), float)
                # str_im = spec.getStringDataArrays()[0]
                # im = np.array([float(s) for s in str_im]).astype(np.float32)
                im_tmp = spec.getFloatDataArrays()[0]
                im = im_tmp.get_data()
                add_df = pd.DataFrame({'native_id': spec.getNativeID(), 'ms_level': spec.getMSLevel(
                ), 'peptide': spec.getMetaValue('peptide'), 'mz': mz, 'rt': rt, 'im': im, 'int': intensity})
                results_df = pd.concat([results_df, add_df])
        return results_df

    @method_timer
    def save_filtered_tsv(self, ms_level=[1], out_file="diapasef_extracted_data.tsv"):
        """
        Save a tsv file of the targeted filtered spectra in tabular format

        Params:
          self: (object) self object that contains filtered data
          ms_level: (list) list of ms level data to write out to file
          out_file: (str) output file to write data to
        """
        results_df = pd.DataFrame()
        for k in tqdm(range(self.filtered.getNrSpectra())):
            spec = self.filtered.getSpectrum(k)
            if spec.getMSLevel() in ms_level:
                mz, intensity = spec.get_peaks()
                if (len(mz) == 0 and len(intensity) == 0):
                    logging.warn(
                        f"MS{spec.getMSLevel()} spectrum native id {spec.getNativeID()} had no m/z or intensity array, skipping this spectrum")
                    continue
                rt = np.full([mz.shape[0]], spec.getRT(), float)
                # str_im = spec.getStringDataArrays()[0]
                # im = np.array([float(s) for s in str_im]).astype(np.float32)
                im_tmp = spec.getFloatDataArrays()[0]
                im = im_tmp.get_data()
                precursor = spec.getPrecursors()[0]
                if self.verbose == 10:
                    logging.debug(
                        f"Adding MS{spec.getMSLevel()} spectrum for peptide: {spec.getMetaValue('peptide')} with native id: {spec.getNativeID()}")
                add_df = pd.DataFrame({'native_id': spec.getNativeID(), 'ms_level': spec.getMSLevel(), 'peptide': spec.getMetaValue(
                    'peptide'), 'precursor_mz': precursor.getMZ(), 'charge': precursor.getCharge(), 'mz': mz, 'rt': rt, 'im': im, 'int': intensity})
                results_df = pd.concat([results_df, add_df])
        logging.info(
            f"Saving filtered target spectra data to tabular tsv format: {out_file}")
        results_df.to_csv(out_file, sep="\t")

    def get_target_ms_level_indices(self):
        """
        Extract spectrum indices for a specific mslevel(s).

        Params:
          self: (object) self object containing meta_data

        Return:
          Return a list of indices with request mslevel(s) to self
        """
        spectra = self.meta_data.getSpectra()
        mslevel_indices = np.array([indice for indice, spec in enumerate(
            spectra) if spec.getMSLevel() in self.mslevel])
        self.mslevel_indices = mslevel_indices

    def get_spectra_rt_list(self):
        """
        Get a list of RT for all the spectra using meta_exp

        Params:
          self: (object) self object containing meta_data

        Return:
          Return a list of RT values for spectra
        """
        meta_rt_list = np.array([meta_spec.getRT()
                                for meta_spec in self.meta_data.getSpectra()])
        self.meta_rt_list = meta_rt_list


class TargeteddiaPASEFExperiment(data_io):
    """
    Class for a targeted DIA-PASEF data extraction
    """

    def __init__(self, mzml_file, peptides, mz_tol=20, rt_window=50, im_window=0.06, mslevel=[1], verbose=0, log_file=None, threads=1):
        """
        Initialize data_access
        """
        super().__init__(mzml_file, mz_tol, rt_window,
                         im_window, mslevel, verbose, log_file)
        self.peptides = peptides
        self.threads = threads

    def set_product(self, mz):
        '''
        Create a product container and set the mz
        '''
        product = po.Product()
        product.setMZ(mz)
        return product

    def get_upper_lower_tol(self, target_mz):
        """
        Get the upper bound and lower bound mz around a target mz given a mz tolerance in ppm

        params:
          target_mz: (float) The target mz to generate upper and lower bound mz coordinates for
          mz_tol: (int) The m/z tolerance toget get upper and lower bounds arround target mz. Must be in ppm.

        Return: 
          (tuple) a tuple of the lower and upper bound for target mz
        """
        # mz_uncertainty = target_mz*(self.mz_tol/(2*1000000))
        mz_uncertainty = target_mz*(self.mz_tol/(1000000))
        target_precursor_mz_upper = target_mz + mz_uncertainty
        target_precursor_mz_lower = target_mz - mz_uncertainty
        return target_precursor_mz_lower, target_precursor_mz_upper

    def get_rt_upper_lower(self, rt_apex):
        """
        Get the upper bound and lower bound of a target RT point for a given RT window, i.e. a window of 50 would be 25 points to either side of the target RT.

        params:
          rt_apex: (float) The target RT point to generate upper and lower bound RT coordinates for
          rt_window: (int) The total window range of RT.

        Return:
          (tuple) a tuple of the lower and upper bound for target RT
        """
        return rt_apex-(self.rt_window/2), rt_apex+(self.rt_window/2)

    def get_im_upper_lower(self, im_apex):
        """
        Get the upper bound and lower bound of a target IM point for a given IM window, i.e. a window of 0.06 would be 0.03 points to either side of the target IM.

        params:
          im_apex: (float) The target IM point to generate upper and lower bound IM coordinates for
          im_window: (int) The total window range of IM.

        Return:
          (tuple) a tuple of the lower and upper bound for target IM
        """
        return im_apex-(self.im_window/2), im_apex+(self.im_window/2)

    def is_mz_in_product_mz_tol_window(self, check_mz, target_product_upper_lower_list):
        """
        Check to see if a target product ion mz is within upper and lower tolerance window

        Params: 
          check_mz: (float) The target product ion mz to check
          target_product_upper_lower_list: (list) A list of tuples that contain an upper and lower bound of a product mz

        Return: (bool) Return a logical value 
        """
        return any([check_mz >= bounds[0] and check_mz <= bounds[1] for bounds in target_product_upper_lower_list])

    # @wrap_non_picklable_objects # seems to cause an error with spectrum. list not having attribute getRT()
    def filter_single_spectrum(self, spec_indice, mslevel, target_peptide_group, rt_start, rt_end, im_start, im_end, target_precursor_mz_lower, target_precursor_mz_upper, target_product_upper_lower_list, verbose=0):
        """
        Filter a single spectrum for a given spectrum indice.

        Params:
          self: (TargeteddiaPASEFExperiment object) an object of self that contains an experiment of OnDiskMSExperiment
          spec_indice: (int) an interger of the spectrum indice to extra a spectrum for
          mslevel: (list) a list of intergers for ms levels to filter for
          target_peptide: (str) target peptide sequence string to set as meta data
          rt_start: (float) start of rt window to filter for
          rt_end: (float) end of rt window to filter for 
          im_start: (float) start of ion mobility window to filter for 
          im_end: (float) end of ion mobility window to filter for 
          target_precursor_mz_lower: (float) the lower acceptable bound of the target precursor within m/z tolerance 
          target_precursor_mz_upper: (float) the upper acceptable bound of the target precursor within m/z tolerance 
          target_product_upper_lower_list: (list) a list of tuples containing lower and upper bound for target product m/z within m/z tolerance
          verbose: (int) an interger specifying verbosity level

        Return:
          If self contains a consumer, filtered spectrum is written to disk, otherwise the filtered spectrum MSSpectrum object is retuned
        """
        spec = self.exp.getSpectrum(spec_indice)
        # TODO: Could remove this if RT check since we already restrict the spectra indices for our given RT range.
        if spec.getRT() >= rt_start and spec.getRT() <= rt_end:
            # Get data arrays
            mz_array = spec.get_peaks()[0]
            int_array = spec.get_peaks()[1]
            im_array = spec.getFloatDataArrays()
            # im_array = im_array[0].get_data() 
            # im_array = [str(im) for im in im_array]
            # im_array = [float(im) for im in im_array]
            check_im_array(im_array[0].get_data() )
            # TODO: sometimes the im_array memory view still gets destroyed or shows different float values?
            im_match_bool = (im_array[0].get_data()  > im_start) & (im_array[0].get_data()  < im_end)
            # TODO: Think about how to vectorize this for-loop. Done.
            if spec.getMSLevel() == 1 and 1 in mslevel:
                mz_match_bool = (mz_array > target_precursor_mz_lower) & (mz_array < target_precursor_mz_upper)
                if verbose == 10 and any(mz_match_bool*im_match_bool):
                    logging.debug(
                        f"Adding MS1 spectrum {spec.getNativeID()} with spectrum indice {spec_indice} filtered for {sum(mz_match_bool*im_match_bool)} spectra between {target_precursor_mz_lower} m/z and {target_precursor_mz_upper} m/z and IM between {im_start} and {im_end}")
            elif spec.getMSLevel() == 2 and 2 in mslevel:
                mz_match_bool = np.array(list(map(self.is_mz_in_product_mz_tol_window, mz_array, itertools.repeat(target_product_upper_lower_list, len(mz_array)) )))
                if verbose == 10 and any(mz_match_bool*im_match_bool):
                    logging.debug(
                        f"Adding MS2 spectrum {spec.getNativeID()} with spectrum indice {spec_indice} filtered for {sum(mz_match_bool*im_match_bool)} spectra between {target_product_upper_lower_list} m/z and IM between {im_start} and {im_end}")
            # Only write out filtered spectra if there is any fitlered spectra to write out
            if any(mz_match_bool*im_match_bool):
                extract_target_indices = np.where(mz_match_bool * im_match_bool)
                filtered_mz = mz_array[extract_target_indices]
                filtered_int = int_array[extract_target_indices]
                filtered_im = im_array[0].get_data() [extract_target_indices]
                # replace peak data with filtered peak data
                spec.set_peaks((filtered_mz, filtered_int))
                # repalce float data arrays with filtered ion mobility data
                fda = po.FloatDataArray()
                filtered_im_np = np.array(filtered_im).astype(np.float32)
                fda.set_data(filtered_im_np)
                fda.setName("Ion Mobility")
                # TODO: There currently is an issue when setting float data array. Getting Float Data Array does not match input
                spec.setFloatDataArrays([fda])
                # Temp solution: Add string data of filtered ion mobility data
                # TODO: Remove temp solution, since it is no longer needed.
                sda = po.StringDataArray()
                sda.setName("String Ion Mobility")
                _ = [sda.push_back(str(im_val)) for im_val in filtered_im]
                spec.setStringDataArrays([sda])
                # Set peptide meta data
                spec.setMetaValue(
                    'peptide', self.peptides[target_peptide_group]['peptide'])

                precursor = po.Precursor()
                precursor.setCharge(self.peptides[target_peptide_group]['charge'])
                precursor.setMZ(
                    self.peptides[target_peptide_group]['precursor_mz'])
                spec.setPrecursors([precursor])

                # Set Prodcut mz values
                spec.setProducts([self.set_product(
                    mz, ) for mz in self.peptides[target_peptide_group]['product_mz']])

                # Write out filtered spectra to file if consumer present. This is more memory efficient
                # TODO: Once paralization works, will need to think about how writing to consumer will be affected
                if self.consumer is not None:
                    self.consumer.consumeSpectrum(spec)
                else:
                    # If you have a lot of filtered spectra to return, it becomes memory heavy.
                    return spec

    @method_timer
    def reduce_spectra(self, output_spectra_file=None):
        """
        Main method for filtering raw mzML diaPASEF data given specific set of coordinates to filter for.

        Params:
          self: (TargeteddiaPASEFExperiment object) an object of self that contains an experiment of OnDiskMSExperiment
          output_spectra_file: (str) output mzML file to save filtered spectra to. Recommended for memory efficiency

        Return:
          if an output file is given, filtered data is written to disk, otherwise a filtered MSExperiment object is returned
        """
        # # Load the data
        # self.load_data()

        # Get Consumer to write filtered spectra to. Memory efficient
        if output_spectra_file is not None:
            logging.debug(
                f"Initializing data consumer to write data to disk - {output_spectra_file}")
            self.get_consumer(output_spectra_file)

        # Get indices for requested ms level spectra
        if self.verbose == 10:
            logging.debug(f"Extracting indices for MS Levels: {self.mslevel}")
        self.get_target_ms_level_indices()

        # Get RT of spectra
        if self.verbose == 10:
            logging.debug(f"Extracting RT values across spectra")
        self.get_spectra_rt_list()
        # logger.info(f"RT Range: {}")

        # TODO: Is this efficient? Or should spectrum access only be extracted during actual spectrum filtering?
        # Main idea behing pre-extracting all desired ms level spectra is to allow for easier paralllization
        # Probably not the best idea to read all spectra into memory. defeats the purpose of on-disk experiment. Will need to find a better way for parallization
        # with code_block_timer(f'Extracting {len(self.mslevel_indices)} of {self.exp.getNrSpectra()} for MS Level {self.mslevel} to spread across {self.threads} threads...'):
        #       spectra_dict = { indice:self.exp.getSpectrum(indice) for indice in self.mslevel_indices.tolist()}

        # Initialize empty MSExperiment container to store filtered data
        if self.consumer is None:
            filtered = po.MSExperiment()
        pbar = tqdm(self.peptides.keys())
        pbar_desc = "INFO: Processing"
        for target_peptide_group in pbar:
            # Update progess bar description
            # pbar_desc = pbar_desc + f"..\n{target_peptide_group}"
            pbar_desc = f"INFO: Processing..{target_peptide_group}"
            pbar.set_description(pbar_desc)

            # Get Coordinates for current peptide
            target_peptide = self.peptides[target_peptide_group]['peptide']
            target_precursor_mz = self.peptides[target_peptide_group]['precursor_mz']
            target_product_mz = self.peptides[target_peptide_group]['product_mz']
            rt_apex = self.peptides[target_peptide_group]['rt_apex']
            im_apex = self.peptides[target_peptide_group]['im_apex']

            # Get tolerance bounds on coordinates
            target_precursor_mz_lower, target_precursor_mz_upper = self.get_upper_lower_tol(
                target_precursor_mz)
            target_product_upper_lower_list = [
                self.get_upper_lower_tol(mz) for mz in target_product_mz]

            rt_start, rt_end = self.get_rt_upper_lower(rt_apex)
            im_start, im_end = self.get_im_upper_lower(im_apex)

            logging.info(
                f"Extracting data for {target_peptide_group} | target precursor: {target_precursor_mz} m/z ({target_precursor_mz_lower} - {target_precursor_mz_upper}) | RT: {rt_apex} sec ({rt_start} - {rt_end}) | IM: {im_apex} 1/Ko ({im_start} - {im_end})")

            # Restrict spectra list further for RT window, to reduce the number of spectra we need to check to perform filtering on
            use_rt_spec_indices = np.where(np.array(
                [spec_rt >= rt_start and spec_rt <= rt_end for spec_rt in self.meta_rt_list]))
            target_spectra_indices = np.intersect1d(
                self.mslevel_indices, use_rt_spec_indices)
            # spectra_list = [spectra_dict[indice] for indice in np.intersect1d(self.mslevel_indices, use_rt_spec_indices)]

            with code_block_timer(f"Filtering {target_spectra_indices.shape[0]} Spectra for {target_peptide_group}..."):
                try:
                    # TODO: Get prallelization to work.
                    # Current Error is with cython __reduce__
                    # File "/home/justincsing/.local/share/r-miniconda/envs/r-reticulate/lib/python3.6/pickle.py", line 496, in save
                    #     rv = reduce(self.proto)
                    #   File "stringsource", line 2, in pyopenms.pyopenms_1.MSSpectrum.__reduce_cython__
                    # TypeError: self.inst cannot be converted to a Python object for pickling
                    # TODO: Currently very memory heavy, maybe save filtered spectrums to disk? Done.
                    filt_spec_list = Parallel(n_jobs=1)(delayed(self.filter_single_spectrum)(spectrum_indice, self.mslevel, target_peptide_group, rt_start, rt_end, im_start, im_end,
                                                                                             target_precursor_mz_lower, target_precursor_mz_upper, target_product_upper_lower_list, self.verbose) for spectrum_indice in target_spectra_indices.tolist())
                except:
                    traceback.print_exc(file=sys.stdout)

            # Add filtered spectra to MSExperiment container
            if self.consumer is None:
                filt_spec_list = list(
                    filter(lambda spectra: spectra is not None, filt_spec_list))
                # Add filtered spectrum to filtered MSExperiment container
                _ = [filtered.addSpectrum(spec) for spec in filt_spec_list]

        if self.consumer is not None:
            # Close consumer to write final part of data out to file
            del self.consumer
        else:
            self.filtered = filtered
