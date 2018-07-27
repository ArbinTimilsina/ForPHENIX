AUTOMAKE_OPTIONS = foreign

INCLUDES = -I$(includedir) -I$(OFFLINE_MAIN)/include -I$(ROOTSYS)/include

lib_LTLIBRARIES = \
  libRun13Jet.la 

noinst_HEADERS = \
  Run13JetLinkDef.h 

include_HEADERS = \
  Run13Jet.h \
  EmcMap.h

libRun13Jet_la_SOURCES = \
  Run13Jet.C \
  Run13Jet_Dict.C

libRun13Jet_la_LIBADD = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -lSubsysReco \
  -lfastjet \
  -lTOAD

nobase_dist_pkgdata_DATA = \
  warnmap.txt

# Rule for generating CINT dictionaries from class headers.
Run13Jet_Dict.C: \
  Run13Jet.h \
  Run13JetLinkDef.h
	rootcint -f $@ -c $(DEFAULT_INCLUDES) $(INCLUDES) $^

clean-local:
	rm -f *Dict*

