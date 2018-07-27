AUTOMAKE_OPTIONS = foreign

INCLUDES = -I. -I$(srcdir) -I$(includedir) -I$(OFFLINE_MAIN)/include -I$(ROOTSYS)/include

lib_LTLIBRARIES = \
  libJetSimAcceptance.la 

noinst_HEADERS = \
  JetSimAcceptanceLinkDef.h 

include_HEADERS = \
  JetSimAcceptance.h 

libJetSimAcceptance_la_SOURCES = \
  JetSimAcceptance.C \
  JetSimAcceptance_Dict.C

libJetSimAcceptance_la_LIBADD = $(OFFLINE_MAIN)/lib -lphool -lphg4hit -lcemc -lPHPythia -L$(ROOTSYS)/lib -lPhysics -lHist -lTree

# Rule for generating CINT dictionaries from class headers.
JetSimAcceptance_Dict.C: \
  JetSimAcceptance.h \
  JetSimAcceptanceLinkDef.h
	rootcint -f $@ -c $(DEFAULT_INCLUDES) $(INCLUDES) $^

clean-local:
	rm -f *Dict*

