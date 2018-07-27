AUTOMAKE_OPTIONS = foreign

INCLUDES = -I$(includedir) -I$(OFFLINE_MAIN)/include -I$(ROOTSYS)/include

lib_LTLIBRARIES = \
  libRun13Jet.la 

noinst_HEADERS = \
  Run13JetLinkDef.h

include_HEADERS = \
  Run13Jet.h \
  EmcMap.h \
  PC3Matching.h \
  TrackQualityPP.h \
  EmcMatching.h \
  EmcTofCorr.h \
  Run13MuCorr.h

libRun13Jet_la_SOURCES = \
  Run13Jet.C \
  Run13Jet_Dict.C \
  PC3Matching.C \
  EmcMatching.C \
  EmcTofCorr.C \
  Run13MuCorr.C

libRun13Jet_la_LIBADD = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -lSubsysReco \
  -lfastjet \
  -luspin \
  -lTOAD

nobase_dist_pkgdata_DATA = \
  warnmap.txt \
  offsetEMC/mean.txt \
  offsetEMC/sigma.txt \
  pointsEMC/Points_Initial_dPhi_Mean.txt \
  pointsEMC/Points_Initial_dPhi_Sigma.txt \
  pointsEMC/Points_Initial_dZ_Mean.txt \
  pointsEMC/Points_Initial_dZ_Sigma.txt \
  pointsEMC/Points_Intermediate_sdPhi_Mean.txt \
  pointsEMC/Points_Intermediate_sdPhi_Sigma.txt \
  pointsEMC/Points_Intermediate_sdZ_Mean.txt \
  pointsEMC/Points_Intermediate_sdZ_Sigma.txt \
  pointsPC3/Points_Initial_dPhi_Mean.txt \
  pointsPC3/Points_Initial_dPhi_Sigma.txt \
  pointsPC3/Points_Initial_dZ_Mean.txt \
  pointsPC3/Points_Initial_dZ_Sigma.txt \
  pointsPC3/Points_Intermediate_sdPhi_Mean.txt \
  pointsPC3/Points_Intermediate_sdPhi_Sigma.txt \
  pointsPC3/Points_Intermediate_sdZ_Mean.txt \
  pointsPC3/Points_Intermediate_sdZ_Sigma.txt \
  emcTOFmap/Run13pp510_EMC_TOF_Correction.root \
  muCorr/lumitree_run13Fin.root

# Rule for generating CINT dictionaries from class headers.
Run13Jet_Dict.C: \
  Run13Jet.h \
  Run13JetLinkDef.h
	rootcint -f $@ -c $(DEFAULT_INCLUDES) $(INCLUDES) $^

clean-local:
	rm -f *Dict*

