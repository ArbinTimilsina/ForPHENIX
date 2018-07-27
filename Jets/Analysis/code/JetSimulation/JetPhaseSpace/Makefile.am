AUTOMAKE_OPTIONS = foreign

INCLUDES = -I. -I$(srcdir) -I$(includedir) -I$(OFFLINE_MAIN)/include -I$(ROOTSYS)/include

lib_LTLIBRARIES = \
  libJetPhaseSpace.la 

noinst_HEADERS = \
  JetPhaseSpaceLinkDef.h 

include_HEADERS = \
  JetPhaseSpace.h

libJetPhaseSpace_la_SOURCES = \
  JetPhaseSpace.C \
  JetPhaseSpace_Dict.C

libJetPhaseSpace_la_LIBADD = $(OFFLINE_MAIN)/lib -lphool -lphg4hit -lcemc -lPHPythia -L$(ROOTSYS)/lib -lPhysics -lHist -lTree

# Rule for generating CINT dictionaries from class headers.
JetPhaseSpace_Dict.C: \
  JetPhaseSpace.h \
  JetPhaseSpaceLinkDef.h
	rootcint -f $@ -c $(DEFAULT_INCLUDES) $(INCLUDES) $^

clean-local:
	rm -f *Dict*

