AUTOMAKE_OPTIONS = foreign

INCLUDES = -I. -I$(srcdir) -I$(includedir) -I$(OFFLINE_MAIN)/include -I$(ROOTSYS)/include

lib_LTLIBRARIES = \
  libJetTriggerPythia.la

noinst_HEADERS = \
  JetTriggerPythiaLinkDef.h

include_HEADERS = \
  JetTriggerPythia.h \
  TAcceptParticle.h

libJetTriggerPythia_la_SOURCES = \
  JetTriggerPythia.C \
  TAcceptParticle.C \
  JetTriggerPythia_Dict.C

libJetTriggerPythia_la_LIBADD = $(OFFLINE_MAIN)/lib -lphool -lPHPythia -L$(ROOTSYS)/lib -lPhysics -lHist -lTree

# Rule for generating CINT dictionaries from class headers.
JetTriggerPythia_Dict.C: \
  JetTriggerPythia.h \
  JetTriggerPythiaLinkDef.h
	rootcint -f $@ -c $(DEFAULT_INCLUDES) $(INCLUDES) $^

clean-local:
	rm -f *Dict*

