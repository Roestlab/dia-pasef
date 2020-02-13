#!/usr/bin/env python
from __future__ import print_function
import pandas as pd
import numpy as np
import sqlite3
import os, sys
from ctypes import *

try:
    if sys.platform[:5] == "win32" or sys.platform[:5] == "win64":
        libname = "timsdata.dll"
    elif sys.platform[:5] == "linux":
        libname = "libtimsdata.so"
    else:
        raise Exception("Unsupported platform.")

    dll = cdll.LoadLibrary(os.path.realpath(libname))
    dll.tims_open.argtypes = [ c_char_p, c_uint32 ]
    dll.tims_open.restype = c_uint64
    dll.tims_close.argtypes = [ c_uint64 ]
    dll.tims_close.restype = None
    dll.tims_get_last_error_string.argtypes = [ c_char_p, c_uint32 ]
    dll.tims_get_last_error_string.restype = c_uint32
    dll.tims_has_recalibrated_state.argtypes = [ c_uint64 ]
    dll.tims_has_recalibrated_state.restype = c_uint32
    dll.tims_read_scans_v2.argtypes = [ c_uint64, c_int64, c_uint32, c_uint32, c_void_p, c_uint32 ]
    dll.tims_read_scans_v2.restype = c_uint32

    convfunc_argtypes = [ c_uint64, c_int64, POINTER(c_double), POINTER(c_double), c_uint32 ]

    dll.tims_index_to_mz.argtypes = convfunc_argtypes
    dll.tims_index_to_mz.restype = c_uint32
    dll.tims_mz_to_index.argtypes = convfunc_argtypes
    dll.tims_mz_to_index.restype = c_uint32

    dll.tims_scannum_to_oneoverk0.argtypes = convfunc_argtypes
    dll.tims_scannum_to_oneoverk0.restype = c_uint32
    dll.tims_oneoverk0_to_scannum.argtypes = convfunc_argtypes
    dll.tims_oneoverk0_to_scannum.restype = c_uint32

    dll.tims_scannum_to_voltage.argtypes = convfunc_argtypes
    dll.tims_scannum_to_voltage.restype = c_uint32
    dll.tims_voltage_to_scannum.argtypes = convfunc_argtypes
    dll.tims_voltage_to_scannum.restype = c_uint32

    def throwLastTimsDataError (dll_handle):
        """Throw last TimsData error string as an exception."""

        len = dll_handle.tims_get_last_error_string(None, 0)
        buf = create_string_buffer(len)
        dll_handle.tims_get_last_error_string(buf, len)
        raise RuntimeError(buf.value)

    print("Found Bruker sdk. Access to the raw data is possible. \n")
    sdk = True
except OSError:
    print("Bruker sdk not found. Some functionalities that need access to raw data will not be available. To activate that functionality place libtimsdata.so (Linux) or timsdata.dll in the src folder. \n")
    sdk = False

class TimsData:


    def get_conversion_func(self):
        q = self.conn.execute("SELECT * FROM TimsCalibration")
        calib = q.fetchone()
        def convert_im(frame_id, im):
            im = np.array(im)
            frame_id = frame_id # Currently the conversion does not depend on the frame_id
            return(1/(calib[8]+calib[9]/(calib[4]+((calib[5]-calib[4])/calib[3])*(im-calib[6]-calib[2]))))
        # Mobility[1/k0] = 1/(c6+c7/(c2+((c3-c2)/c1)*(scanno-c4-c0)))
        return convert_im

    if sdk:

        def __init__ (self, analysis_directory, use_recalibrated_state=False):
            if sys.version_info.major == 2:
                if not isinstance(analysis_directory, unicode):
                    raise ValueError("analysis_directory must be a Unicode string.")
            if sys.version_info.major == 3:
                if not isinstance(analysis_directory, str):
                    raise ValueError("analysis_directory must be a string.")

            self.dll = dll

            self.handle = self.dll.tims_open(
                analysis_directory.encode('utf-8'),
                1 if use_recalibrated_state else 0 )
            if self.handle == 0:
                throwLastTimsDataError(self.dll)

            self.conn = sqlite3.connect(os.path.join(analysis_directory, "analysis.tdf"))

            q = self.conn.execute("SELECT MsMsType FROM Frames")

            COLUMN = 0
            frames = q.fetchall()
            ms2_ids = [row[COLUMN] for row in frames]
            self.ms2_id = max(ms2_ids)

            self.initial_frame_buffer_size = 128 # may grow in readScans()
            self.dll = dll

        def __del__ (self):
            if hasattr(self, 'handle'):
                self.dll.tims_close(self.handle)

        def __callConversionFunc (self, frame_id, input_data, func):
            if type(input_data) is np.ndarray and input_data.dtype == np.float64:
                # already "native" format understood by DLL -> avoid extra copy
                in_array = input_data
            else:
                # convert data to format understood by DLL:
                in_array = np.array(input_data, dtype=np.float64)

            cnt = len(in_array)
            out = np.empty(shape=cnt, dtype=np.float64)
            success = func(self.handle, frame_id,
                           in_array.ctypes.data_as(POINTER(c_double)),
                           out.ctypes.data_as(POINTER(c_double)),
                           cnt)
            if success == 0:
                throwLastTimsDataError(self.dll)
            return out

        def indexToMz (self, frame_id, mzs):
            return self.__callConversionFunc(frame_id, mzs, self.dll.tims_index_to_mz)

        def mzToIndex (self, frame_id, mzs):
            return self.__callConversionFunc(frame_id, mzs, self.dll.tims_mz_to_index)

        def scanNumToOneOverK0 (self, frame_id, mzs):
            return self.__callConversionFunc(frame_id, mzs, self.dll.tims_scannum_to_oneoverk0)

        def oneOverK0ToScanNum (self, frame_id, mzs):
            return self.__callConversionFunc(frame_id, mzs, self.dll.tims_oneoverk0_to_scannum)

        def scanNumToVoltage (self, frame_id, mzs):
            return self.__callConversionFunc(frame_id, mzs, self.dll.tims_scannum_to_voltage)

        def voltageToScanNum (self, frame_id, mzs):
            return self.__callConversionFunc(frame_id, mzs, self.dll.tims_voltage_to_scannum)


        # Output: list of tuples (indices, intensities)
        def readScans (self, frame_id, scan_begin, scan_end):
            # buffer-growing loop
            while True:
                cnt = int(self.initial_frame_buffer_size) # necessary cast to run with python 3.5
                buf = np.empty(shape=cnt, dtype=np.uint32)
                len = 4 * cnt

                required_len = self.dll.tims_read_scans_v2(self.handle, frame_id, scan_begin, scan_end,
                                                            buf.ctypes.data_as(POINTER(c_uint32)),
                                                            len)
                if required_len == 0:
                    throwLastTimsDataError(self.dll)

                if required_len > len:
                    if required_len > 16777216:
                        # arbitrary limit for now...
                        raise RuntimeError("Maximum expected frame size exceeded.")
                    self.initial_frame_buffer_size = required_len / 4 + 1 # grow buffer
                else:
                    break

            result = []
            d = scan_end - scan_begin
            for i in range(scan_begin, scan_end):
                npeaks = buf[i-scan_begin]
                indices     = buf[d : d+npeaks]
                d += npeaks
                intensities = buf[d : d+npeaks]
                d += npeaks
                result.append((indices,intensities))

            return result

        def readSpectrum (self, frame_id, scan_begin, scan_end):

            scans = self.readScans(frame_id, scan_begin, scan_end)
            # Summarize on a grid
            allind = []
            allint = []
            for scan in scans:
                indices = np.array(scan[0])
                if len(indices) > 0:
                    # summed_intensities = np.zeros(indices.max() + 1)
                    intens = scan[1]
                    allind = np.concatenate((allind, indices))
                    allint = np.concatenate((allint, intens))
                    # summed_intensities[indices] += intens
            allmz = self.indexToMz(frame_id, allind)
            return((allmz, allint))

    else:

        def __init__ (self, analysis_directory, use_recalibrated_state=False):
            if sys.version_info.major == 2:
                if not isinstance(analysis_directory, unicode):
                    raise ValueError("analysis_directory must be a Unicode string.")
            if sys.version_info.major == 3:
                if not isinstance(analysis_directory, str):
                    raise ValueError("analysis_directory must be a string.")

            self.conn = sqlite3.connect(os.path.join(analysis_directory, "analysis.tdf"))

        def scanNumToOneOverK0 (self, frame_id, mzs):
            convert_scan_num = self.get_conversion_func()
            return(convert_scan_num(frame_id, mzs))

class DiaPasefData(TimsData):

    def get_cycle_length (self):
        firstcycle = True
        initial_ms1 = True
        cycle_length = 0
        row = 1
        while (firstcycle):
            mstype = self.conn.execute("SELECT MsMsType FROM Frames LIMIT {0},1".format(row)).fetchone()
            row += 1
            if mstype[0] == self.ms2_id:
                cycle_length += 1
                initial_ms1 = False
            elif(initial_ms1 != True):
                firstcycle = False
        return cycle_length

    def get_windows_per_frame(self):
        if self.ms2_id == 8:
            q = self.conn.execute("SELECT * FROM Frames INNER JOIN PasefFrameMsMsInfo ON Frames.Id=PasefFrameMsMsInfo.Frame WHERE MsMsType={0} AND Frame = (SELECT MIN(Frame) FROM PasefFrameMsMsInfo as fr)".format(self.ms2_id))
        elif self.ms2_id == 9:
            q = self.conn.execute("SELECT * FROM Frames INNER JOIN DiaFrameMsMsInfo ON Frames.Id = DiaFrameMsMsInfo.Frame INNER JOIN DiaFrameMsMsWindows ON DiaFrameMsMsInfo.WindowGroup = DiaFrameMsMsWindows.WindowGroup WHERE MsMsType={0} AND Frame = (SELECT MIN(Frame) FROM DiaFrameMsMsInfo as fr)".format(self.ms2_id))
        windows = q.fetchall()
        return len(windows)

    def get_windows(self):
        """Extracts the window scheme from the first cycle of a tims file"""
        #pasef_ms2_id = 8 # diaPASEF ms2 scans are denoted by 8 instead of 2
        if self.ms2_id == 8:
            cycle_length = self.get_cycle_length()
            wpf = self.get_windows_per_frame()
            q = self.conn.execute("SELECT * FROM Frames INNER JOIN PasefFrameMsMsInfo ON Frames.Id=PasefFrameMsMsInfo.Frame WHERE MsMsType={0} LIMIT {1}" .format(self.ms2_id, cycle_length*wpf))
        elif self.ms2_id == 9:
            cycle_length = self.get_cycle_length()
            wpf = self.get_windows_per_frame()
            q = self.conn.execute("SELECT * FROM Frames INNER JOIN DiaFrameMsMsInfo ON Frames.Id = DiaFrameMsMsInfo.Frame INNER JOIN DiaFrameMsMsWindows ON DiaFrameMsMsInfo.WindowGroup = DiaFrameMsMSWindows.WindowGroup WHERE MsMsType={0} LIMIT {1}".format(self.ms2_id, cycle_length*wpf))
        frames = q.fetchall()
        colnames = [description[0] for description in q.description]
        resframe = pd.DataFrame(data = frames, columns = colnames)
        # Currently the frame does not matter for the conversion, thus we can set it to 1
        resframe["IMstart"] = self.scanNumToOneOverK0(1, resframe.ScanNumBegin)
        resframe["IMend"] = self.scanNumToOneOverK0(1, resframe.ScanNumEnd)
        return resframe
