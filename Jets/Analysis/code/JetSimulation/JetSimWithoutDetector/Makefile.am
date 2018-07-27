AUTOMAKE_OPTIONS = foreign

INCLUDES = -I. -I$(srcdir) -I$(includedir) -I$(OFFLINE_MAIN)/include -I$(ROOTSYS)/include

lib_LTLIBRARIES = \
  libJetSimWithoutDetector.la 

noinst_HEADERS = \
  JetSimWithoutDetectorLinkDef.h 

include_HEADERS = \
  JetSimWithoutDetector.h \
  TAcceptParticle.h

libJetSimWithoutDetector_la_SOURCES = \
  JetSimWithoutDetector.C \
  TAcceptParticle.C \
  JetSimWithoutDetector_Dict.C

libJetSimWithoutDetector_la_LIBADD = $(OFFLINE_MAIN)/lib -lphool -lphg4hit -lcemc -lPHPythia -L$(ROOTSYS)/lib -lPhysics -lHist -lTree -lHepMC

# Rule for generating CINT dictionaries from class headers.
JetSimWithoutDetector_Dict.C: \
  JetSimWithoutDetector.h \
  JetSimWithoutDetectorLinkDef.h
	rootcint -f $@ -c $(DEFAULT_INCLUDES) $(INCLUDES) $^

clean-local:
	rm -f *Dict*

