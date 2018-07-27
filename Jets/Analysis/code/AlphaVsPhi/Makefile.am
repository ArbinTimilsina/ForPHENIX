 AUTOMAKE_OPTIONS = foreign

 INCLUDES = -I$(includedir) -I$(OFFLINE_MAIN)/include -I$(ROOTSYS)/include 

 lib_LTLIBRARIES = libAlphaVsPhi.la

 noinst_HEADERS = \
   AlphaVsPhiLinkDef.h
 
 include_HEADERS = \
   AlphaVsPhi.h

 AM_LDFLAGS = \
   -L$(libdir) \
   -L$(OFFLINE_MAIN)/lib

 libAlphaVsPhi_la_LIBADD = \
   -L$(libdir) \
   -lCNT

 libAlphaVsPhi_la_SOURCES = \
   AlphaVsPhi.C \
   AlphaVsPhi_Dict.C
 
 BUILT_SOURCES = \
   testexternals.C

 noinst_PROGRAMS = \
   testexternals

 testexternals_LDADD = \
   libAlphaVsPhi.la

 testexternals.C:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

 AlphaVsPhi_Dict.C: AlphaVsPhi.h AlphaVsPhiLinkDef.h
	rootcint -f $@ -c $(DEFAULT_INCLUDES) $(INCLUDES) $^
