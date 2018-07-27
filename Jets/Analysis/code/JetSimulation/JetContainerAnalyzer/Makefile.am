AUTOMAKE_OPTIONS = foreign

INCLUDES = -I. -I$(srcdir) -I$(includedir) -I$(OFFLINE_MAIN)/include -I$(ROOTSYS)/include

lib_LTLIBRARIES = \
  libJetContainerAnalyzer.la 

noinst_HEADERS = \
  JetContainerAnalyzerLinkDef.h 

include_HEADERS = \
  JetContainerAnalyzer.h \
  TrackQualityCuAu.h \
  TrackQualityPP.h \
  TAcceptParticle.h

libJetContainerAnalyzer_la_SOURCES = \
  JetContainerAnalyzer.C \
  TAcceptParticle.C \
  JetContainerAnalyzer_Dict.C

libJetContainerAnalyzer_la_LIBADD = $(OFFLINE_MAIN)/lib -lphool -lPHPythia -L$(ROOTSYS)/lib -lPhysics -lHist -lTree

# Rule for generating CINT dictionaries from class headers.
JetContainerAnalyzer_Dict.C: \
  JetContainerAnalyzer.h \
  JetContainerAnalyzerLinkDef.h
	rootcint -f $@ -c $(DEFAULT_INCLUDES) $(INCLUDES) $^

clean-local:
	rm -f *Dict*

