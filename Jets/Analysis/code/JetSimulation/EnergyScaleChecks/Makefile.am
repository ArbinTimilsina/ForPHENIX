AUTOMAKE_OPTIONS = foreign

INCLUDES = -I. -I$(srcdir) -I$(includedir) -I$(OFFLINE_MAIN)/include -I$(ROOTSYS)/include

lib_LTLIBRARIES = \
  libEnergyScaleChecks.la 

noinst_HEADERS = \
  EnergyScaleChecksLinkDef.h 

include_HEADERS = \
  EnergyScaleChecks.h 

libEnergyScaleChecks_la_SOURCES = \
  EnergyScaleChecks.C \
  EnergyScaleChecks_Dict.C

libEnergyScaleChecks_la_LIBADD = $(OFFLINE_MAIN)/lib -lphool -lPHPythia -L$(ROOTSYS)/lib -lPhysics -lHist -lTree

# Rule for generating CINT dictionaries from class headers.
EnergyScaleChecks_Dict.C: \
  EnergyScaleChecks.h \
  EnergyScaleChecksLinkDef.h
	rootcint -f $@ -c $(DEFAULT_INCLUDES) $(INCLUDES) $^

clean-local:
	rm -f *Dict*

