
# HeLa file analysis

First, we split the file

    time python convert_single.py ../data/hela/ /tmp/hela_all_frames.np.mzML # takes ca 2.5 h

Then we split into individual mzML files (per SWATH):

    OpenSwathFileSplitter -in  /tmp/hela_all_frames.np.mzML -outputDirectory /tmp/hela_all_split/ # takes ca 1.5 h

Then we extract chromatograms:

    for f in /tmp/hela_all_split/openswath_tmpfile_*.mzML; 
      do echo $f;
      bn=`basename $f .mzML`;
      echo $bn;
      OpenSwathChromatogramExtractor -in  $f -tr /media/data/tmp/CiRT_ALL.TraML -out /tmp/hela_chromout_n_$bn.mzML -mz_window 25 -ppm -is_swath ;
    done
    FileMerger -in /tmp/hela_chromout_n* -out /tmp/hela_allchroms.mzML

Then we search for matching peaks for alignment:

    OpenSwathRTNormalizer  -tr /media/data/tmp/CiRT_ALL.TraML -in /tmp/hela_allchroms.mzML -out /tmp/hela.trafoXML -min_coverage 0.1 -algorithm:TransitionGroupPicker:PeakPickerMRM:peak_width 5
    Rscript plot_trafo.R /tmp/hela.trafoXML  out.pdf

    for f in /tmp/hela_all_split/openswath_tmpfile_*.mzML; 
      do echo $f;
      bn=`basename $f .mzML`;
      echo $bn;
      OpenSwathChromatogramExtractor -in  $f -tr /media/data/tmp/CiRT_ALL.TraML -out /tmp/hela_chromout_n_$bn.mzML -mz_window 25 -ppm -is_swath -rt_norm /tmp/hela.trafoXML -rt_window 1200 
    done
    FileMerger -in /tmp/hela_chromout_align* -out /tmp/hela_allchroms_align.mzML




