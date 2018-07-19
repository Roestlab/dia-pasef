import pyopenms
import sys
import numpy as np
inp = "/tmp/test.mzML"

inp = sys.argv[1]
out = sys.argv[2]
merge_scans = int(sys.argv[3])
drop_scans = int(sys.argv[4])

class MergeConsumer():
    """
        Merging consumer that merges MS2 spectra with the same precursor. The
        number of consecutive spectra to be merged is a user parameter.

        This class merges m/z and intensity coordinates as well as the first
        FloatDataArray.

        Does not implement chromatogram consuming.
    """

    def __init__(self, consumer, merge_nr, drop_nr):
        assert drop_nr <= merge_nr, "Drop number needs to be smaller than merge number"
        self._internal_consumer = consumer
        self._spectrum_storage = {}
        self._merge_nr = merge_nr
        self._drop_nr = drop_nr # how many scans to drop

    def __del__(self):
        """
        Cleanup: write all stored spectra to disk
        """
        for k, v in self._spectrum_storage.items():
            if v:
                # for i, kk in enumerate(v):
                #     print i, kk
                merge_spec = self._mergeSpectra(v)
                self._internal_consumer.consumeSpectrum(merge_spec)

    def setExpectedSize(self, a, b):
        pass

    def setExperimentalSettings(self, s):
        pass

    def consumeChromatogram(self, s):
        self._internal_consumer.consumeChromatogram(c)

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
                tmp_data = tmp[self._drop_nr:]
                # print ("will merge", len(tmp_data), "droping", self._drop_nr )
                merge_spec = self._mergeSpectra(tmp_data)
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

        import pyopenms
        import numpy as np
        merge_spec = pyopenms.MSSpectrum(tmp[0])
        fda = tmp[0].getFloatDataArrays()[0]
        fda.clear()
        allmz = []
        allint = []
        # print ("_mergeSpectr", len(tmp))
        for q in tmp:
            m, i = q.get_peaks()
            allmz.append(m)
            allint.append(i)

            # Sanity checks, precursors of merged spectra and float arrays need to match
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

# print("will start here")
consumer = pyopenms.PlainMSDataWritingConsumer(out)
consumer = MergeConsumer(consumer, merge_scans, drop_scans)

exp = pyopenms.MSExperiment()
pyopenms.MzMLFile().transform(inp, consumer)


