AUTOMAKE_OPTIONS = foreign

INCLUDES = -I. -I$(srcdir) -I$(includedir) -I$(OFFLINE_MAIN)/include -I$(ROOTSYS)/include

lib_LTLIBRARIES = \
  libRun13JetSim.la 

noinst_HEADERS = \
  Run13JetSimLinkDef.h 

include_HEADERS = \
  Run13JetSim.h

libRun13JetSim_la_SOURCES = \
  Run13JetSim.C \
  Run13JetSim_Dict.C

libRun13JetSim_la_LIBADD = $(OFFLINE_MAIN)/lib -lphool -lPHPythia -L$(ROOTSYS)/lib -lPhysics -lHist -lTree

# Rule for generating CINT dictionaries from class headers.
Run13JetSim_Dict.C: \
  Run13JetSim.h \
  Run13JetSimLinkDef.h
	rootcint -f $@ -c $(DEFAULT_INCLUDES) $(INCLUDES) $^

clean-local:
	rm -f *Dict*

