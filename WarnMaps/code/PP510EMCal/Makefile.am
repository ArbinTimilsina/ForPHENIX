AUTOMAKE_OPTIONS = foreign

INCLUDES = -I$(includedir) -I$(OFFLINE_MAIN)/include -I$(ROOTSYS)/include

lib_LTLIBRARIES = \
  libEMCalMapPP510.la 

noinst_HEADERS = \
  EMCalMapPP510LinkDef.h 

include_HEADERS = \
  EMCalMapPP510.h

libEMCalMapPP510_la_SOURCES = \
  EMCalMapPP510.C \
  EMCalMapPP510_Dict.C

libEMCalMapPP510_la_LIBADD = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -lSubsysReco

# Rule for generating CINT dictionaries from class headers.
EMCalMapPP510_Dict.C: \
  EMCalMapPP510.h \
  EMCalMapPP510LinkDef.h
	rootcint -f $@ -c $(DEFAULT_INCLUDES) $(INCLUDES) $^

clean-local:
	rm -f *Dict*

