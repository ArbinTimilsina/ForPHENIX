AUTOMAKE_OPTIONS = foreign

INCLUDES = -I. -I$(srcdir) -I$(includedir) -I$(OFFLINE_MAIN)/include -I$(ROOTSYS)/include

lib_LTLIBRARIES = \
  libJetConstituentPerformance.la 

noinst_HEADERS = \
  JetConstituentPerformanceLinkDef.h 

include_HEADERS = \
  JetConstituentPerformance.h \
  TrackQualityPP.h \
  TAcceptParticle.h

libJetConstituentPerformance_la_SOURCES = \
  JetConstituentPerformance.C \
  TAcceptParticle.C \
  JetConstituentPerformance_Dict.C

libJetConstituentPerformance_la_LIBADD = $(OFFLINE_MAIN)/lib -lphool -lphg4hit -lcemc -lPHPythia -L$(ROOTSYS)/lib -lPhysics -lHist -lTree

# Rule for generating CINT dictionaries from class headers.
JetConstituentPerformance_Dict.C: \
  JetConstituentPerformance.h \
  JetConstituentPerformanceLinkDef.h
	rootcint -f $@ -c $(DEFAULT_INCLUDES) $(INCLUDES) $^

clean-local:
	rm -f *Dict*