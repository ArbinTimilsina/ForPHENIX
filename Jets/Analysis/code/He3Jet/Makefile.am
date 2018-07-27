AUTOMAKE_OPTIONS = foreign

INCLUDES = -I$(includedir) -I$(OFFLINE_MAIN)/include -I$(ROOTSYS)/include

lib_LTLIBRARIES = \
  libHe3Jet.la 

noinst_HEADERS = \
  He3JetLinkDef.h 

include_HEADERS = \
  He3Jet.h \
  EmcMap.h

libHe3Jet_la_SOURCES = \
  He3Jet.C \
  He3Jet_Dict.C

libHe3Jet_la_LIBADD = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -lSubsysReco \
  -lfastjet \
  -lTOAD

nobase_dist_pkgdata_DATA = \
  warnmap.txt

# Rule for generating CINT dictionaries from class headers.
He3Jet_Dict.C: \
  He3Jet.h \
  He3JetLinkDef.h
	rootcint -f $@ -c $(DEFAULT_INCLUDES) $(INCLUDES) $^

clean-local:
	rm -f *Dict*

