AUTOMAKE_OPTIONS = foreign

INCLUDES = -I$(includedir) -I$(OFFLINE_MAIN)/include -I$(ROOTSYS)/include

lib_LTLIBRARIES = \
  libJetAnalyzer.la 

noinst_HEADERS = \
  JetAnalyzerLinkDef.h 

include_HEADERS = \
  JetAnalyzer.h \
  EmcMap.h \
  TrackQualityCuAu.h \
  TrackQualityPP.h \
  Pc3Matching.h \
  EmcMatching.h

libJetAnalyzer_la_SOURCES = \
  JetAnalyzer.C \
  Pc3Matching.C \
  EmcMatching.C \
  JetAnalyzer_Dict.C

libJetAnalyzer_la_LIBADD = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -lSubsysReco \
  -lfastjet \
  -lTOAD

nobase_dist_pkgdata_DATA = \
  warnmapCuAu.txt \
  warnmapPP.txt \
  offset/mean.txt \
  offset/sigma.txt \
  parametersEMC/Parameters_Initial_dPhi_Mean.txt \
  parametersEMC/Parameters_Initial_dPhi_Sigma.txt \
  parametersEMC/Parameters_Initial_dZ_Mean.txt \
  parametersEMC/Parameters_Initial_dZ_Sigma.txt \
  parametersEMC/Parameters_Final_dPhi_Mean.txt \
  parametersEMC/Parameters_Final_dPhi_Sigma.txt \
  parametersEMC/Parameters_Final_dZ_Mean.txt \
  parametersEMC/Parameters_Final_dZ_Sigma.txt \
  parametersPC3/CuAu200_pc3_par_C0_I.txt \
  parametersPC3/CuAu200_pc3_par_C0_II.txt \
  parametersPC3/CuAu200_pc3_par_C1_I.txt \
  parametersPC3/CuAu200_pc3_par_C1_II.txt \
  parametersPC3/CuAu200_pc3_par_C2_I.txt \
  parametersPC3/CuAu200_pc3_par_C2_II.txt \
  parametersPC3/CuAu200_pc3_par_C3_I.txt \
  parametersPC3/CuAu200_pc3_par_C3_II.txt \
  parametersPC3/CuAu200_pc3_par_C4_I.txt \
  parametersPC3/CuAu200_pc3_par_C4_II.txt \
  parametersPC3/CuAu200_pc3_par_C5_I.txt \
  parametersPC3/CuAu200_pc3_par_C5_II.txt \
  parametersPC3/CuAu200_pc3_par_C6_I.txt \
  parametersPC3/CuAu200_pc3_par_C6_II.txt \
  parametersPC3/CuAu200_pc3_par_C7_I.txt \
  parametersPC3/CuAu200_pc3_par_C7_II.txt \
  parametersPC3/CuAu200_pc3_par_C8_I.txt \
  parametersPC3/CuAu200_pc3_par_C8_II.txt \
  parametersPC3/CuAu200_pc3_par_C9_I.txt \
  parametersPC3/CuAu200_pc3_par_C9_II.txt \
  parametersPC3/CuAu200_pc3_slopes_C0_I.txt \
  parametersPC3/CuAu200_pc3_slopes_C0_II.txt \
  parametersPC3/CuAu200_pc3_slopes_C1_I.txt \
  parametersPC3/CuAu200_pc3_slopes_C1_II.txt \
  parametersPC3/CuAu200_pc3_slopes_C2_I.txt \
  parametersPC3/CuAu200_pc3_slopes_C2_II.txt \
  parametersPC3/CuAu200_pc3_slopes_C3_I.txt \
  parametersPC3/CuAu200_pc3_slopes_C3_II.txt \
  parametersPC3/CuAu200_pc3_slopes_C4_I.txt \
  parametersPC3/CuAu200_pc3_slopes_C4_II.txt \
  parametersPC3/CuAu200_pc3_slopes_C5_I.txt \
  parametersPC3/CuAu200_pc3_slopes_C5_II.txt \
  parametersPC3/CuAu200_pc3_slopes_C6_I.txt \
  parametersPC3/CuAu200_pc3_slopes_C6_II.txt \
  parametersPC3/CuAu200_pc3_slopes_C7_I.txt \
  parametersPC3/CuAu200_pc3_slopes_C7_II.txt \
  parametersPC3/CuAu200_pc3_slopes_C8_I.txt \
  parametersPC3/CuAu200_pc3_slopes_C8_II.txt \
  parametersPC3/CuAu200_pc3_slopes_C9_I.txt \
  parametersPC3/CuAu200_pc3_slopes_C9_II.txt

# Rule for generating CINT dictionaries from class headers.
JetAnalyzer_Dict.C: \
  JetAnalyzer.h \
  JetAnalyzerLinkDef.h
	rootcint -f $@ -c $(DEFAULT_INCLUDES) $(INCLUDES) $^

clean-local:
	rm -f *Dict*

