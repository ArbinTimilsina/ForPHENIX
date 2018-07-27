AUTOMAKE_OPTIONS = foreign

INCLUDES = -I. -I$(srcdir) -I$(includedir) -I$(OFFLINE_MAIN)/include -I$(ROOTSYS)/include

lib_LTLIBRARIES = \
  libJetHadCorrection.la 

noinst_HEADERS = \
  JetHadCorrectionLinkDef.h

include_HEADERS = \
  JetHadCorrection.h

libJetHadCorrection_la_SOURCES = \
  JetHadCorrection.C \
  JetHadCorrection_Dict.C

libJetHadCorrection_la_LIBADD = $(OFFLINE_MAIN)/lib -lphool -lPHPythia -L$(ROOTSYS)/lib -lPhysics -lHist -lTree

# Rule for generating CINT dictionaries from class headers.
JetHadCorrection_Dict.C: \
  JetHadCorrection.h \
  JetHadCorrectionLinkDef.h
	rootcint -f $@ -c $(DEFAULT_INCLUDES) $(INCLUDES) $^

clean-local:
	rm -f *Dict*

