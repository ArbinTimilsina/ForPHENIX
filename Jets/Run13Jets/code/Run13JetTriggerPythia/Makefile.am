AUTOMAKE_OPTIONS = foreign

INCLUDES = -I. -I$(srcdir) -I$(includedir) -I$(OFFLINE_MAIN)/include -I$(ROOTSYS)/include

lib_LTLIBRARIES = \
  libRun13JetTriggerPythia.la

noinst_HEADERS = \
  Run13JetTriggerPythiaLinkDef.h

include_HEADERS = \
  Run13JetTriggerPythia.h

libRun13JetTriggerPythia_la_SOURCES = \
  Run13JetTriggerPythia.C \
  Run13JetTriggerPythia_Dict.C

libRun13JetTriggerPythia_la_LIBADD = $(OFFLINE_MAIN)/lib -lphool -lPHPythia -L$(ROOTSYS)/lib -lPhysics -lHist -lTree

# Rule for generating CINT dictionaries from class headers.
Run13JetTriggerPythia_Dict.C: \
  Run13JetTriggerPythia.h \
  Run13JetTriggerPythiaLinkDef.h
	rootcint -f $@ -c $(DEFAULT_INCLUDES) $(INCLUDES) $^

clean-local:
	rm -f *Dict*

