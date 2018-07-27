AUTOMAKE_OPTIONS = foreign

INCLUDES = -I. -I$(srcdir) -I$(includedir) -I$(OFFLINE_MAIN)/include -I$(ROOTSYS)/include

lib_LTLIBRARIES = \
  libJetNloCorrection.la 

noinst_HEADERS = \
  JetNloCorrectionLinkDef.h

include_HEADERS = \
  JetNloCorrection.h \
  TAcceptParticle.h

libJetNloCorrection_la_SOURCES = \
  JetNloCorrection.C \
  TAcceptParticle.C \
  JetNloCorrection_Dict.C

libJetNloCorrection_la_LIBADD = $(OFFLINE_MAIN)/lib -lphool -lphg4hit -lcemc -lPHPythia -L$(ROOTSYS)/lib -lPhysics -lHist -lTree -lHepMC

# Rule for generating CINT dictionaries from class headers.
JetNloCorrection_Dict.C: \
  JetNloCorrection.h \
  JetNloCorrectionLinkDef.h
	rootcint -f $@ -c $(DEFAULT_INCLUDES) $(INCLUDES) $^

clean-local:
	rm -f *Dict*

