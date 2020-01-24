#!/usr/bin/env python
from __future__ import print_function
import pyopenms
import numpy as np


class MergeConsumer(object):
    """
        Merging consumer that merges MS2 spectra with the same precursor. The
        number of consecutive spectra to be merged is a user parameter.

        This class merges m/z and intensity coordinates as well as the first
        FloatDataArray.

        Does not implement chromatogram consuming.
    """

    def __init__(self, consumer, merge_nr, check_precursor = True, reference_spectrum = 0):
        self._internal_consumer = consumer
        self._spectrum_storage = {}
        self._merge_nr = merge_nr

        self._check_precursor = check_precursor
        self._reference_spectrum = reference_spectrum

    def __del__(self):
        """
        Cleanup: write all stored spectra to disk
        """
        for k, v in self._spectrum_storage.items():
            if v:
                merge_spec = self._mergeSpectra(v)
                self._internal_consumer.consumeSpectrum(merge_spec)

    def consumeSpectrum(self, s):
        """
            consume individual spectra:
                - write MS1 spectra directly to disk
                - collect MS2 spectra and write afer merging n spectra together
        """
        if s.getMSLevel() == 1:
            self._internal_consumer.consumeSpectrum(s)
        else:

            mz = int(s.getPrecursors()[0].getMZ()*10)
            tmp = self._spectrum_storage.get(mz, [])
            tmp.append(s)

            if len(tmp) >= self._merge_nr:
                merge_spec = self._mergeSpectra(tmp)
                self._internal_consumer.consumeSpectrum(merge_spec)
                tmp = []

            self._spectrum_storage[mz] = tmp

    def _mergeSpectra(self, tmp):
        """
        Perform spectral merging
         - we pick one spectrum (the first one) as our reference and add
           all other data from the other spectra to it
         - we check that merged spectra have equal precursor m/z and same float array
        """

        ref_spec = self._reference_spectrum
        if ref_spec >= len(tmp): 
            ref_spec = 0
        merge_spec = pyopenms.MSSpectrum(tmp[ref_spec])
        fda = tmp[ref_spec].getFloatDataArrays()[0]
        
        fda.clear()
        allmz = []
        allint = []
        for q in tmp:
            m, i = q.get_peaks()
            allmz.append(m)
            allint.append(i)

            # Sanity checks, precursors of merged spectra and float arrays need to match
            if self._check_precursor:
                assert q.getPrecursors()[0].getMZ() - merge_spec.getPrecursors()[0].getMZ() < 1e-5
                assert len(q.getFloatDataArrays()) == len(merge_spec.getFloatDataArrays())
                assert q.getFloatDataArrays()[0].getName() == merge_spec.getFloatDataArrays()[0].getName()

            # TODO this is not very efficient, fix in pyOpenMS!
            for d in q.getFloatDataArrays()[0]:
                fda.push_back(d)

        mz = np.concatenate(allmz)
        intens = np.concatenate(allint)

        # create merge spec
        merge_spec.set_peaks( (mz, intens) )
        merge_spec.setFloatDataArrays([fda])
        merge_spec.sortByPosition()
        return merge_spec

class TenzerMergeConsumer(MergeConsumer):
    """
        Merging consumer that merges any set of N consecutive MS2 spectra. The
        number of consecutive spectra to be merged is a user parameter.

        This class merges m/z and intensity coordinates as well as the first
        FloatDataArray.

        Does not implement chromatogram consuming.
    """

    def __init__(self, consumer, merge_nr, total_scans, check_precursor = False, reference_spectrum = 0):
        super(TenzerMergeConsumer, self).__init__(consumer, merge_nr, check_precursor, reference_spectrum)
        self._cnt = 0
        self._total_nr = total_scans

    def consumeSpectrum(self, s):
        """
            consume individual spectra:
                - write MS1 spectra directly to disk
                - collect MS2 spectra and write afer merging n spectra together
        """
        if s.getMSLevel() == 1:
            self._internal_consumer.consumeSpectrum(s)
        else:

            # Iterate through all internal storage lists to which we need to
            # append the current spectrum
            mz = int(s.getPrecursors()[0].getMZ())
            for k in range(self._cnt - self._merge_nr + 1, self._cnt + 1, 1):
                # Skip those lists that are at the edge
                if k < 0: continue
                if k > self._total_nr - self._merge_nr: continue

                tmp = self._spectrum_storage.get(k, [])
                tmp.append(s)

                if len(tmp) >= self._merge_nr:
                    merge_spec = self._mergeSpectra(tmp)
                    self._internal_consumer.consumeSpectrum(merge_spec)
                    mz = int(merge_spec.getPrecursors()[0].getMZ())
                    tmp = []

                self._spectrum_storage[k] = tmp

            self._cnt += 1
            self._cnt %= self._total_nr

