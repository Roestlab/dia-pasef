#!/usr/bin/env python
from __future__ import print_function

class SplittingConsumer():
    """
        Splitting consumer that deals with splitting up MS2 spectra from the same m/z window.

        This class will distribute incoming spectra into N different downstream
        consumers based on the incoming order and precursor m/z. This can be
        useful when dealing with overlapping windows where spectra need to be
        distributed into alternating files (e.g. the first spectrum with
        precursor 400-425 will go into the first consumer, the second spectrum
        with precursor 400-425 will go into the second consumer, the next one
        into the first consumer again etc.).

        Does not implement chromatogram consuming.
    """

    def __init__(self, consumers):
        self._internal_consumers = consumers
        self._nr_splits = len(consumers)
        self._curr_consumers = {}
        self._cnt = 0

    def __del__(self):
        """
        Cleanup: write all stored spectra to disk
        """
        pass

    def consumeSpectrum(self, s):
        """
            consume individual spectra:
                - write MS1 spectra to all consumers
                - identify current  MS2 spectra and write afer merging n spectra together
        """
        if s.getMSLevel() == 1:
            for c in self._internal_consumers:
                c.consumeSpectrum(s)
        else:

            mz = int(s.getPrecursors()[0].getMZ()*10)
            tmp = self._curr_consumers.get(mz, 0)

            c = self._internal_consumers[tmp]
            c.consumeSpectrum(s)

            self._cnt = self._cnt + 1

            # identify the next consumer to be used
            next_consumer = (tmp+1) % self._nr_splits
            self._curr_consumers[mz] = next_consumer

