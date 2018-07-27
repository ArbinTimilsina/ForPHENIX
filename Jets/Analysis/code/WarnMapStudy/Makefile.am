AUTOMAKE_OPTIONS = foreign

INCLUDES = -I$(includedir) -I$(OFFLINE_MAIN)/include -I$(ROOTSYS)/include

lib_LTLIBRARIES = \
  libWarnMapStudy.la 

noinst_HEADERS = \
  WarnMapStudyLinkDef.h 

include_HEADERS = \
  WarnMapStudy.h \
  EmcMap.h

libWarnMapStudy_la_SOURCES = \
  WarnMapStudy.C \
  EmcMap.C \
  WarnMapStudy_Dict.C

libWarnMapStudy_la_LIBADD = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -lSubsysReco \
  -lTOAD

nobase_dist_pkgdata_DATA = \
  warnmap.txt

# Rule for generating CINT dictionaries from class headers.
WarnMapStudy_Dict.C: \
  WarnMapStudy.h \
  WarnMapStudyLinkDef.h
	rootcint -f $@ -c $(DEFAULT_INCLUDES) $(INCLUDES) $^

clean-local:
	rm -f *Dict*

