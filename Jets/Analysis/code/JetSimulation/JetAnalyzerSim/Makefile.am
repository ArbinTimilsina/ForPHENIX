AUTOMAKE_OPTIONS = foreign

INCLUDES = -I. -I$(srcdir) -I$(includedir) -I$(OFFLINE_MAIN)/include -I$(ROOTSYS)/include

lib_LTLIBRARIES = \
  libJetAnalyzerSim.la 

noinst_HEADERS = \
  JetAnalyzerSimLinkDef.h 

include_HEADERS = \
  JetAnalyzerSim.h

libJetAnalyzerSim_la_SOURCES = \
  JetAnalyzerSim.C \
  JetAnalyzerSim_Dict.C

libJetAnalyzerSim_la_LIBADD = $(OFFLINE_MAIN)/lib -lphool -lPHPythia -L$(ROOTSYS)/lib -lPhysics -lHist -lTree

# Rule for generating CINT dictionaries from class headers.
JetAnalyzerSim_Dict.C: \
  JetAnalyzerSim.h \
  JetAnalyzerSimLinkDef.h
	rootcint -f $@ -c $(DEFAULT_INCLUDES) $(INCLUDES) $^

clean-local:
	rm -f *Dict*

