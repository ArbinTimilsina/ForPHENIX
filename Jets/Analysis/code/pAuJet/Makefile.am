AUTOMAKE_OPTIONS = foreign

INCLUDES = -I$(includedir) -I$(OFFLINE_MAIN)/include -I$(ROOTSYS)/include

lib_LTLIBRARIES = \
  libpAuJet.la 

noinst_HEADERS = \
  pAuJetLinkDef.h 

include_HEADERS = \
  pAuJet.h \
  EmcMap.h

libpAuJet_la_SOURCES = \
  pAuJet.C \
  pAuJet_Dict.C

libpAuJet_la_LIBADD = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -lSubsysReco \
  -lfastjet \
  -lTOAD

nobase_dist_pkgdata_DATA = \
  warnmap.txt

# Rule for generating CINT dictionaries from class headers.
pAuJet_Dict.C: \
  pAuJet.h \
  pAuJetLinkDef.h
	rootcint -f $@ -c $(DEFAULT_INCLUDES) $(INCLUDES) $^

clean-local:
	rm -f *Dict*

