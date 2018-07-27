AUTOMAKE_OPTIONS = foreign

INCLUDES = -I. -I$(srcdir) -I$(includedir) -I$(OFFLINE_MAIN)/include -I$(ROOTSYS)/include

lib_LTLIBRARIES = \
  libJetTriggerData.la

noinst_HEADERS = \
  JetTriggerDataLinkDef.h

include_HEADERS = \
  JetTriggerData.h

libJetTriggerData_la_SOURCES = \
  JetTriggerData.C \
  JetTriggerData_Dict.C

libJetTriggerData_la_LIBADD = $(OFFLINE_MAIN)/lib -lphool -L$(ROOTSYS)/lib -lPhysics -lHist -lTree

# Rule for generating CINT dictionaries from class headers.
JetTriggerData_Dict.C: \
  JetTriggerData.h \
  JetTriggerDataLinkDef.h
	rootcint -f $@ -c $(DEFAULT_INCLUDES) $(INCLUDES) $^

clean-local:
	rm -f *Dict*

