AUTOMAKE_OPTIONS = foreign

INCLUDES = -I$(includedir) -I$(OFFLINE_MAIN)/include -I$(ROOTSYS)/include

lib_LTLIBRARIES = \
  libEMCalMap.la 

noinst_HEADERS = \
  EMCalMapLinkDef.h 

include_HEADERS = \
  EMCalMap.h

libEMCalMap_la_SOURCES = \
  EMCalMap.C \
  EMCalMap_Dict.C

libEMCalMap_la_LIBADD = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -lSubsysReco

# Rule for generating CINT dictionaries from class headers.
EMCalMap_Dict.C: \
  EMCalMap.h \
  EMCalMapLinkDef.h
	rootcint -f $@ -c $(DEFAULT_INCLUDES) $(INCLUDES) $^

clean-local:
	rm -f *Dict*

