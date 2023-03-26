import os
from pathlib import Path
from typing import Optional
import pandas as pd
import numpy as np
import sqlite3
import zlib
import traceback
import logging
import multiprocessing 
from tqdm import tqdm
from datetime import datetime

def get_upper_lower_tol(target_mz_list: list, mz_tol: int):
    """
    Get the upper bound and lower bound mz around a target mz given a mz tolerance in ppm

    params:
      target_mz_list: (list) The list of target mz to generate upper and lower bound mz coordinates for
      mz_tol: (int) The m/z tolerance toget get upper and lower bounds arround target mz. Must be in ppm.

    Return: 
      (tuple) a tuple of the lower and upper bound for target mz
    """
    target_mz_list = np.asarray(target_mz_list)
    mz_uncertainty = target_mz_list * mz_tol / 2.0 * 1.0e-6
    target_precursor_mz_upper = target_mz_list + mz_uncertainty
    target_precursor_mz_lower = target_mz_list - mz_uncertainty
    return target_precursor_mz_lower, target_precursor_mz_upper

def create_database(db_filename: str) -> None:
    """
    Create a new SQLite database with the specified filename and schema.

    Parameters:
        db_filename (str): The name of the new database file to be created.

    Returns:
        None
    """
    # Create a connection to the database
    conn = sqlite3.connect(db_filename)

    # Create the tables
    conn.execute('''CREATE TABLE DATA(SPECTRUM_ID INT,CHROMATOGRAM_ID INT,MOBILOGRAM_ID INT, RAWEXTRACT_ID INT,COMPRESSION INT,DATA_TYPE INT,DATA BLOB NOT NULL)''')
    conn.execute('''CREATE TABLE SPECTRUM(ID INT PRIMARY KEY NOT NULL,RUN_ID INT,MSLEVEL INT NULL,RETENTION_TIME REAL NULL,SCAN_POLARITY INT NULL,NATIVE_ID TEXT NOT NULL)''')
    conn.execute('''CREATE TABLE RUN(ID INT PRIMARY KEY NOT NULL,FILENAME TEXT NOT NULL, NATIVE_ID TEXT NOT NULL)''')
    conn.execute('''CREATE TABLE RUN_EXTRA(RUN_ID INT,DATA BLOB NOT NULL)''')
    conn.execute('''CREATE TABLE CHROMATOGRAM(ID INT PRIMARY KEY NOT NULL,RUN_ID INT,NATIVE_ID TEXT NOT NULL)''')
    conn.execute('''CREATE TABLE MOBILOGRAM(ID INT PRIMARY KEY NOT NULL,RUN_ID INT,NATIVE_ID TEXT NOT NULL)''')
    conn.execute('''CREATE TABLE RAW_EXTRACTION(ID INT PRIMARY KEY NOT NULL,RUN_ID INT,NATIVE_ID TEXT NOT NULL)''')
    conn.execute('''CREATE TABLE PRODUCT(SPECTRUM_ID INT,CHROMATOGRAM_ID INT,MOBILOGRAM_ID INT, RAWEXTRACT_ID INT,CHARGE INT NULL,ISOLATION_TARGET REAL NULL,ISOLATION_LOWER REAL NULL,ISOLATION_UPPER REAL NULL)''')
    conn.execute('''CREATE TABLE PRECURSOR(SPECTRUM_ID INT,CHROMATOGRAM_ID INT,MOBILOGRAM_ID INT, RAWEXTRACT_ID INT,CHARGE INT NULL,PEPTIDE_SEQUENCE TEXT NULL,DRIFT_TIME REAL NULL,ACTIVATION_METHOD INT NULL,ACTIVATION_ENERGY REAL NULL,ISOLATION_TARGET REAL NULL,ISOLATION_LOWER REAL NULL,ISOLATION_UPPER REAL NULL)''')

    # Create the indexes
    conn.execute('''CREATE INDEX data_chr_idx ON DATA(CHROMATOGRAM_ID)''')
    conn.execute('''CREATE INDEX data_sp_idx ON DATA(SPECTRUM_ID)''')
    conn.execute('''CREATE INDEX data_mob_idx ON DATA(MOBILOGRAM_ID)''')
    conn.execute('''CREATE INDEX data_raw_idx ON DATA(RAWEXTRACT_ID)''')
    conn.execute('''CREATE INDEX spec_rt_idx ON SPECTRUM(RETENTION_TIME)''')
    conn.execute('''CREATE INDEX spec_mslevel_idx ON SPECTRUM(MSLEVEL)''')
    conn.execute('''CREATE INDEX spec_run_idx ON SPECTRUM(RUN_ID)''')
    conn.execute('''CREATE INDEX run_extra_idx ON RUN_EXTRA(RUN_ID)''')
    conn.execute('''CREATE INDEX chrom_run_idx ON CHROMATOGRAM(RUN_ID)''')
    conn.execute('''CREATE INDEX mobi_run_idx ON MOBILOGRAM(RUN_ID)''')
    conn.execute('''CREATE INDEX raw_run_idx ON RAW_EXTRACTION(RUN_ID)''')
    conn.execute('''CREATE INDEX product_chr_idx ON DATA(CHROMATOGRAM_ID)''')
    conn.execute('''CREATE INDEX product_sp_idx ON DATA(SPECTRUM_ID)''')
    conn.execute('''CREATE INDEX product_mob_idx ON DATA(MOBILOGRAM_ID)''')
    conn.execute('''CREATE INDEX product_raw_idx ON DATA(RAWEXTRACT_ID)''')
    conn.execute('''CREATE INDEX precursor_chr_idx ON DATA(CHROMATOGRAM_ID)''')
    conn.execute('''CREATE INDEX precursor_sp_idx ON DATA(SPECTRUM_ID)''')
    conn.execute('''CREATE INDEX precursor_mob_idx ON DATA(MOBILOGRAM_ID)''')
    conn.execute('''CREATE INDEX precursor_raw_idx ON DATA(RAWEXTRACT_ID)''')

    # Commit the changes and close the connection
    conn.commit()
    conn.close()

def insert_meta_mapping_data(id: int, pep_coord: dict, prec_info: pd.DataFrame, mslevel: int = 1) -> str:
    """
    Returns SQL queries to insert precursor, product, spectrum, chromatogram and mobilogram data into the database.

    Args:
        id (int): Unique identifier for the data to be inserted.
        pep_coord (dict): A dict containing precursor and product ion information for a peptide.
        prec_info (pd.DataFrame): A pandas DataFrame containing precursor ion information.
        mslevel (int, optional): The MS level of the data. Defaults to 1.

    Returns:
        str: SQL queries to insert precursor, product, spectrum, chromatogram and mobilogram data into the database.
    """
    # Create a list to store the SQL queries
    sql_queries = []

    # Insert Precursor data
    prec_meta = prec_info.loc[prec_info.ms_level.isin([mslevel])]
    precursor_values = (id, id, id, id, pep_coord['charge'], pep_coord['peptide'], pep_coord['assay_im'], 1, 0.0, pep_coord['precursor_mz'], 0.0, 0.0)
    insert_precursor_query = f"INSERT INTO PRECURSOR (SPECTRUM_ID, CHROMATOGRAM_ID, MOBILOGRAM_ID, RAWEXTRACT_ID, CHARGE, PEPTIDE_SEQUENCE, DRIFT_TIME, ACTIVATION_METHOD, ACTIVATION_ENERGY, ISOLATION_TARGET, ISOLATION_LOWER, ISOLATION_UPPER) VALUES ({precursor_values[0]}, {precursor_values[1]}, {precursor_values[2]}, {precursor_values[3]}, {precursor_values[4]}, '{precursor_values[5]}', {precursor_values[6]}, {precursor_values[7]}, {precursor_values[8]}, {precursor_values[9]}, {precursor_values[10]}, {precursor_values[11]});"
    sql_queries.append(insert_precursor_query)

    ## Insert Product data
    if mslevel == 1:
        product_values = (id, id, id, id, 0, 0.0, 0.0, 0.0)
    elif mslevel == 2:
        product_values = (id, id, id, id, pep_coord['product_charge'][np.where(pep_coord['product_mz']==prec_info.mz_annotation.values[0])[0][0]], prec_info.mz_annotation.values[0], 0.0, 0.0)
    insert_product_query = "INSERT INTO PRODUCT (SPECTRUM_ID, CHROMATOGRAM_ID, MOBILOGRAM_ID, RAWEXTRACT_ID, CHARGE, ISOLATION_TARGET, ISOLATION_LOWER, ISOLATION_UPPER) VALUES ({}, {}, {}, {}, {}, {}, {}, {});".format(*product_values)
    sql_queries.append(insert_product_query)

    # Insert Spectrum data
    if mslevel == 1:
        spectrum_values = (id, pep_coord['run_id'], 1, pep_coord['rt_apex'], 0, str(pep_coord['precursor_id']) + "_" + "Precursor")
    elif mslevel == 2:
        spectrum_values = (id, pep_coord['run_id'], 2, pep_coord['rt_apex'], 0, pep_coord['transition_id'][np.where(pep_coord['product_mz']==prec_info.mz_annotation.values[0])[0][0]])
    insert_spectrum_query = f"INSERT INTO SPECTRUM (ID, RUN_ID, MSLEVEL, RETENTION_TIME, SCAN_POLARITY, NATIVE_ID) VALUES ({spectrum_values[0]}, {spectrum_values[1]}, {spectrum_values[2]}, {spectrum_values[3]}, {spectrum_values[4]}, '{spectrum_values[5]}');"
    sql_queries.append(insert_spectrum_query)

    ## Insert Chromatogram data
    if mslevel == 1:
        chromatogram_values = (id, pep_coord['run_id'], str(pep_coord['precursor_id']) + "_" + "Precursor")
    elif mslevel == 2:
        chromatogram_values = (id, pep_coord['run_id'], pep_coord['transition_id'][np.where(pep_coord['product_mz']==prec_info.mz_annotation.values[0])[0][0]])
    insert_chromatogram_query = f"INSERT INTO CHROMATOGRAM (ID, RUN_ID, NATIVE_ID) VALUES ({chromatogram_values[0]}, {chromatogram_values[1]}, '{chromatogram_values[2]}');"
    sql_queries.append(insert_chromatogram_query)

    ## Insert Mobilogram data
    if mslevel == 1:
        mobilogram_values = (id, pep_coord['run_id'], str(pep_coord['precursor_id']) + "_" + "Precursor")
    elif mslevel == 2:
        mobilogram_values = (id, pep_coord['run_id'], pep_coord['transition_id'][np.where(pep_coord['product_mz']==prec_info.mz_annotation.values[0])[0][0]])
    insert_mobilogram_query = f"INSERT INTO MOBILOGRAM (ID, RUN_ID, NATIVE_ID) VALUES ({mobilogram_values[0]}, {mobilogram_values[1]}, '{mobilogram_values[2]}');"
    sql_queries.append(insert_mobilogram_query)

    ## Insert Raw Extraction data
    if mslevel == 1:
        mobilogram_values = (id, pep_coord['run_id'], str(pep_coord['precursor_id']) + "_" + "Precursor")
    elif mslevel == 2:
        mobilogram_values = (id, pep_coord['run_id'], pep_coord['transition_id'][np.where(pep_coord['product_mz']==prec_info.mz_annotation.values[0])[0][0]])
    insert_mobilogram_query = f"INSERT INTO RAW_EXTRACTION (ID, RUN_ID, NATIVE_ID) VALUES ({mobilogram_values[0]}, {mobilogram_values[1]}, '{mobilogram_values[2]}');"
    sql_queries.append(insert_mobilogram_query)

    # Combine all the queries into a transaction string
    transaction_query = ""
    for query in sql_queries:
        transaction_query += query + "\n"
    return transaction_query

def insert_compressed_data(data: pd.DataFrame, id: int, compression_level: int, data_type: int, data_id: int) -> str:
    """
     Builds a SQL query for compressing and inserting the given data into the 'DATA' table.

    - compression is one of 0 = no, 1 = zlib, 2 = np-linear, 3 = np-slof, 4 = np-pic, 5 = np-linear + zlib, 6 = np-slof + zlib, 7 = np-pic + zlib
    - data_type is one of 0 = mz, 1 = int, 2 = rt, 3 = im, 4 = raw extracted 4D data
    - data contains the raw (blob) data for a single data array

    Parameters:
        data (pd.DataFrame): The data to be compressed and inserted.
        id (int): The ID of the spectrum or chromatogram to which the data belongs.
        compression_level (int): The compression level to be used for compressing the data.
        data_type (int): The type of data (e.g., mz, rt, im, etc.).
        data_id (int): The type of ID (spectrum, chromatogram, or mobilogram) that the data belongs to.

    Returns:
        str: The SQL query for inserting the compressed data into the database.
    """
    if compression_level==1:
        compressed_data = zlib.compress(data.tobytes())
    else:
        compressed_data = data.tobytes()
    # Decompressing back
    # np.frombuffer(zlib.decompress(compressed_data), dtype=np.float64)
    data_to_insert = (id if data_id == 0 else 'NULL', id if data_id == 2 else 'NULL', id if data_id == 3 else 'NULL', id if data_id == 4 else 'NULL', compression_level, data_type, compressed_data.hex())
    # query = "INSERT INTO DATA (SPECTRUM_ID, CHROMATOGRAM_ID, MOBILOGRAM_ID, COMPRESSION, DATA_TYPE, DATA) VALUES ({},{},{},{},{},{});".format(*data_to_insert)
    query = f"INSERT INTO DATA (SPECTRUM_ID, CHROMATOGRAM_ID, MOBILOGRAM_ID, RAWEXTRACT_ID, COMPRESSION, DATA_TYPE, DATA) VALUES ({data_to_insert[0]},{data_to_insert[1]},{data_to_insert[2]},{data_to_insert[3]},{data_to_insert[4]},{data_to_insert[5]}, x'{data_to_insert[6]}');"

    return query

def process_data(data_type: str, group: pd.DataFrame, id: int) -> None:
    """Process data based on the specified data type and insert compressed data into the database.

    Args:
        data_type (str): The type of data to process ('mz', 'rt', or 'im').
        group (pd.DataFrame): The group of data to process.
        id (int): The ID for the data.
        conn (sqlite3.Connection): The connection to the database.

    Raises:
        ValueError: If an invalid data type is specified.

    Returns:
        str: The SQL query for inserting the compressed data into the database.
    """
    # Create a list to store the SQL queries
    sql_queries = []
    if data_type == "mz":
        sorted_group = group.sort_values('mz')
        sorted_group = sorted_group.round({'mz': 4}).groupby(["ms_level", "precursor_id", "peptide",  "precursor_mz",  "charge", "mz"]).agg({'int': np.sum}).reset_index()
        #### Insert mz data
        sql_queries.append(insert_compressed_data(sorted_group['mz'].values, id, 1, 0, 0))
        #### Insert Int data
        sql_queries.append(insert_compressed_data(sorted_group['int'].values, id, 1, 1, 0))
    elif data_type == "rt":
        sorted_group = group.sort_values('rt')
        sorted_group = sorted_group.groupby(["ms_level", "precursor_id", "peptide",  "precursor_mz",  "charge", "rt"]).agg({'int': np.sum}).reset_index()
        #### Insert rt data
        sql_queries.append(insert_compressed_data(sorted_group['rt'].values, id, 1, 2, 2))
        #### Insert Int data
        sql_queries.append(insert_compressed_data(sorted_group['int'].values, id, 1, 1, 2))
    elif data_type == "im":
        sorted_group = group.sort_values('im')
        sorted_group = sorted_group.groupby(["ms_level", "precursor_id", "peptide",  "precursor_mz",  "charge", "im"]).agg({'int': np.sum}).reset_index()
        #### Insert im data
        sql_queries.append(insert_compressed_data(sorted_group['im'].values, id, 1, 3, 3))
        #### Insert Int data
        sql_queries.append(insert_compressed_data(sorted_group['int'].values, id, 1, 1, 3))
    else:
        raise ValueError("Invalid data type")
    transaction_query = ""
    for query in sql_queries:
        transaction_query += query + "\n"
    return transaction_query

def process_precursor(precursor_str: str, coords_data: dict, dt: pd.DataFrame, shared_id: multiprocessing.Manager, lock: multiprocessing.Lock, mz_tol: int = 25) -> None:
    """
    Process precursor data and insert into database.

    Args:
        precursor_str: A string representing the precursor ID.
        coords_data: A dictionary containing precursor coordinates data.
        dt: A pandas DataFrame containing precursor data.
        shared_id: A multiprocessing manager's shared value object used to assign unique IDs for each precursor.
        lock: A multiprocessing manager's lock object to ensure thread safety while accessing shared resources.

    Returns:
        None
    """
    dt_filt = dt.loc[dt['precursor_id'] == precursor_str]
    pep_coord = coords_data[precursor_str]

    dt_filt.loc[:, 'mz_annotation'] = dt_filt['mz']
    if dt_filt.ms_level.isin([1]).any():
        dt_filt.loc[dt_filt.ms_level.isin([1]), ['mz_annotation']] = dt_filt['precursor_mz']

    if dt_filt.ms_level.isin([2]).any():
        try:
            # Extract the 'product_mz' list from the 'pep_coord' dict and convert it to a NumPy array
            prod_mz = np.array(pep_coord['product_mz'])

            # Get the upper and lower tolerances for the product m/z values based on the specified 'mz_tol' value
            prod_lower_mz, prod_upper_mz = get_upper_lower_tol(pep_coord['product_mz'], mz_tol)

            # Create a boolean mask that is True for all rows in 'dt_filt' where the difference between the 'mz' value and 
            # each value in 'prod_mz' is within the range of the upper and lower tolerances, and where the 'ms_level' value 
            # is not equal to 1. 
            mask = (np.abs(dt_filt.mz.values[:, np.newaxis] - prod_mz[:, np.newaxis].transpose()) <= prod_upper_mz - prod_lower_mz) & (dt_filt.ms_level.values[:, np.newaxis]!=1)

            # Select all rows in 'dt_filt' where any element in the corresponding row of the 'mask' array is True, and update 
            # the 'mz_annotation' column for these rows to be the corresponding value from the 'prod_mz' array.
            # dt_filt.loc[mask.any(axis=1), ['mz_annotation']] = prod_mz[np.where(mask)[1]]
            # Note: Use argmax instead of where to find the first True, because there may be multiple True values in the same row, (when there are product ions with the same or similar m/z value)
            dt_filt.loc[mask.any(axis=1), ['mz_annotation']] = prod_mz[np.argmax(mask[mask.any(axis=1)], axis=1)]

        except ValueError as e:
            logging.error(f"Error: {e} for precursor {precursor_str}")
            traceback.print_exc()
            raise e
    
    sql_query_list = []
    try:
        prec_info = dt_filt[['ms_level', 'peptide', 'charge', 'mz_annotation']].drop_duplicates().sort_values('ms_level')
        if any(prec_info.ms_level.isin([1])):
            with lock:
                id = shared_id.value
                shared_id.value += 1
            ## Insert Meta Data Tables
            prec_meta = prec_info.loc[prec_info.ms_level.isin([1])]
            sql_query_list.append(insert_meta_mapping_data(id, pep_coord, prec_meta, 1))
            

            ## Insert Data
            datams1 = dt_filt.loc[dt_filt.ms_level.isin([1])]
            for mz, group in datams1.groupby('mz_annotation'):
                ### Spectrum data
                sql_query_list.append(process_data("mz", group, id))
                ### Chromatogram data
                sql_query_list.append(process_data("rt", group, id))
                ### Mobilogram data
                sql_query_list.append(process_data("im", group, id))
            
            ## Insert Raw Extracted Data
            raw_data = datams1[['mz_annotation', 'mz', 'rt', 'im', 'int']].to_numpy()
            ### Compressing and Decompressing Data
            # compressed_data = zlib.compress(raw_data.tobytes())
            # # decompress the data with zlib
            # decompressed_data = zlib.decompress(compressed_data)
            # # convert the decompressed data to a numpy array
            # data = np.frombuffer(decompressed_data, dtype=np.float64).reshape((-1, 5))
            ### Raw data
            sql_query_list.append(insert_compressed_data(raw_data, id, 1, 4, 4))

            # increment id
            with lock:
                id = shared_id.value
                shared_id.value += 1
    except ValueError as e:
        logging.error(f"Error: {e} for precursor {precursor_str}")
        traceback.print_exc()
        raise e

    try:
        if any(prec_info.ms_level.isin([2])):
            datams2 = dt_filt.loc[dt_filt.ms_level.isin([2])]
            for mz, prod_group in datams2.groupby('mz_annotation'):
                ## Insert Meta Data Tables
                prod_info = prod_group[['ms_level', 'peptide', 'charge', 'mz_annotation']].drop_duplicates().sort_values('ms_level')
                sql_query_list.append(insert_meta_mapping_data(id, pep_coord, prod_info, 2))

                ## Insert Data
                ### Spectrum data
                sql_query_list.append(process_data("mz", prod_group, id))
                ### Chromatogram data
                sql_query_list.append(process_data("rt", prod_group, id))
                ### Mobilogram data
                sql_query_list.append(process_data("im", prod_group, id))

                ## Insert Raw Extracted Data
                raw_data = prod_group[['mz_annotation', 'mz', 'rt', 'im', 'int']].to_numpy()
                sql_query_list.append(insert_compressed_data(raw_data, id, 1, 4, 4))

                # increment id  
                with lock:  
                    id = shared_id.value
                    shared_id.value += 1
    except ValueError as e:
        logging.error(f"Error: {e} for precursor {precursor_str}")
        traceback.print_exc()
        raise e

    transaction_query = ""
    for query in sql_query_list:
        transaction_query += query + "\n"
    return transaction_query

def export_sqmass(file: str, coordsfile: str, db_filename: Optional[str] = None, mz_tol: int = 25):
    """
    Exports precursor mass data to an SQLite database file.

    Parameters:
    -----------
    file : str
        Path to the input file containing precursor mass data in TSV format.
    coordsfile : str
        Path to the input file containing precursor coordinates data in pickle format.
    db_filename : str, optional
        Path to the SQLite database file. If not specified, a default file name will be created based on the input file name.

    Raises:
    -------
    ValueError:
        If the input file or coordinates file is not found.

    Returns:
    --------
    None
    """
    # Check that the input files exist
    file_path = Path(file)
    coords_path = Path(coordsfile)
    if not file_path.is_file():
        raise ValueError(f"File not found: {file}")
    if not coords_path.is_file():
        raise ValueError(f"File not found: {coordsfile}")
    
    if db_filename is None:
        db_filename = os.path.splitext(file)[0] + ".sqMass"

    # Load data
    dt = pd.read_csv(file, sep="\t")
    dt['precursor_id'] = dt['peptide'].apply(str) + '_' + dt['charge'].apply(str)

    # Load precursor coords dictionary
    coords_data = pd.read_pickle(coordsfile)

    # Create database if it doesn't exist
    if os.path.exists(db_filename):
        logging.error(f"Info: Database file {db_filename} already exists.  Skipping creation.") 
    else:
        logging.info(f"Info: Creating database file {db_filename}")
        create_database(db_filename)

    # Create a connection to the database
    conn = sqlite3.connect(db_filename)

    ## Insert Run data
    conn.execute(f"""INSERT INTO RUN (ID, FILENAME, NATIVE_ID) VALUES ({coords_data[list(coords_data.keys())[0]]['run_id']}, '{file}', '{file}')""")
    conn.commit()
    conn.close()

    # precursors
    precursors = np.unique(dt.precursor_id)

    # Set the number of processes to use
    num_processes = 10

    # Disable the SettingWithCopyWarning
    pd.options.mode.chained_assignment = None  # default='warn'

    start_date = datetime.now().strftime("%Y-%m-%d %H:%M:%S") # get the current date and time
    logging.info(f"Started generating queries {start_date}.")
    with multiprocessing.Manager() as manager:
        shared_id = manager.Value('i', 0)
        lock = manager.Lock()
        # Create a pool of processes to write the data to the database
        with multiprocessing.Pool(processes=num_processes) as pool:
            args = [(precursor_str, coords_data, dt, shared_id, lock, mz_tol) for precursor_str in precursors]
            results = list(tqdm(pool.starmap(process_precursor, args), total=len(args)))
    end_date = datetime.now().strftime("%Y-%m-%d %H:%M:%S") # get the current date and time
    logging.info(f"Finished generatign queries {end_date}.")

    # Re-enable the SettingWithCopyWarning
    pd.options.mode.chained_assignment = 'warn'

    # connect to database
    conn = sqlite3.connect(db_filename)
    cursor = conn.cursor()

    transaction_query = "BEGIN TRANSACTION;"
    for query in results:
        if query is not None:
            transaction_query += query + "\n"
    transaction_query += "COMMIT;"

    start_date = datetime.now().strftime("%Y-%m-%d %H:%M:%S") # get the current date and time

    # execute queries using executemany()
    logging.info(f"Started writing data to database at {start_date}.")
    cursor.executescript(transaction_query)
    conn.commit()
    end_date = datetime.now().strftime("%Y-%m-%d %H:%M:%S") # get the current date and time
    logging.info(f"Finished writing data to database at {end_date}.")

    # close connection
    conn.close()