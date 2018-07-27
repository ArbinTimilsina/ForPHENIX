AUTOMAKE_OPTIONS = foreign

INCLUDES = -I$(includedir) -I$(OFFLINE_MAIN)/include -I$(ROOTSYS)/include

lib_LTLIBRARIES = \
  libRunQA.la 

noinst_HEADERS = \
  RunQALinkDef.h 

include_HEADERS = \
  RunQA.h

libRunQA_la_SOURCES = \
  RunQA.C \
  RunQA_Dict.C

libRunQA_la_LIBADD = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -lSubsysReco

# Rule for generating CINT dictionaries from class headers.
RunQA_Dict.C: \
  RunQA.h \
  RunQALinkDef.h
	rootcint -f $@ -c $(DEFAULT_INCLUDES) $(INCLUDES) $^

clean-local:
	rm -f *Dict*

