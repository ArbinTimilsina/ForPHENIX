# -*- mode: makefile; -*-

lib_LTLIBRARIES = libjetrec.la

#####################################################################

libjetrec_la_SOURCES = background.cc recbaseline.cc recext.cc \
	recfilter.cc recfilterc.cc recfilteroptimize.cc \
	recfiltersseiir.cc recfilterssemap.cc recfiltersseopt.cc \
	solve.cc

CXXFLAGS = @CXXFLAGS@ -I@srcdir@ -I@srcdir@/../libjetbase \
	-I@srcdir@/../libjetevent -DHAVE_SSE -DHAVE_SSE2

libjetrec_la_LDFLAGS = @ROOT_LDFLAGS@
libjetrec_la_LIBADD = @ROOT_LIBS@

noinst_HEADERS =
nobase_include_HEADERS = jetrec/background.h jetrec/cuda.h jetrec/jet.h \
	jetrec/rec.h jetrec/series.h jetrec/solve.h jetrec/util.h
pkgdata_DATA =

%_dict.cc:	jetrec/%.h %_linkdef.h
	$(ROOTCINT) -f $@ -c $(DEFS) $(DEFAULT_INCLUDES) \
	$(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) $(AM_CXXFLAGS) \
	$(CXXFLAGS) $^

recfiltersseiir.o:	recfiltersseiir.cc
	g++ $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) \
	$(AM_CPPFLAGS) $(CPPFLAGS) $(AM_CXXFLAGS) $(CXXFLAGS) \
	-O2 -c $< -o $@

recfilterssemap.o:	recfilterssemap.cc
	g++ $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) \
	$(AM_CPPFLAGS) $(CPPFLAGS) $(AM_CXXFLAGS) $(CXXFLAGS) \
	-O2 -c $< -o $@

recfiltersseopt.o:	recfiltersseopt.cc
	g++ $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) \
	$(AM_CPPFLAGS) $(CPPFLAGS) $(AM_CXXFLAGS) $(CXXFLAGS) \
	-O2 -c $< -o $@

if LIBTOOL_1_X
recfiltersseiir.lo:	recfiltersseiir.cc
	g++ $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) \
	$(AM_CPPFLAGS) $(CPPFLAGS) $(AM_CXXFLAGS) $(CXXFLAGS) \
	-O2 -c $< -fPIC -DPIC -o $@

recfilterssemap.lo:	recfilterssemap.cc
	g++ $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) \
	$(AM_CPPFLAGS) $(CPPFLAGS) $(AM_CXXFLAGS) $(CXXFLAGS) \
	-O2 -c $< -fPIC -DPIC -o $@

recfiltersseopt.lo:	recfiltersseopt.cc
	g++ $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) \
	$(AM_CPPFLAGS) $(CPPFLAGS) $(AM_CXXFLAGS) $(CXXFLAGS) \
	-O2 -c $< -fPIC -DPIC -o $@
else
recfiltersseiir.lo:     recfiltersseiir.cc
	/bin/sh ../libtool --tag=CXX --mode=compile g++ \
	$(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) \
        $(AM_CPPFLAGS) $(CPPFLAGS) $(AM_CXXFLAGS) $(CXXFLAGS) \
        -O2 -c $< -fPIC -DPIC -o $@

recfilterssemap.lo:     recfilterssemap.cc
	/bin/sh ../libtool --tag=CXX --mode=compile g++ \
	$(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) \
	$(AM_CPPFLAGS) $(CPPFLAGS) $(AM_CXXFLAGS) $(CXXFLAGS) \
	-O2 -c $< -fPIC -DPIC -o $@

recfiltersseopt.lo:     recfiltersseopt.cc
	/bin/sh ../libtool --tag=CXX --mode=compile g++ \
	$(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) \
	$(AM_CPPFLAGS) $(CPPFLAGS) $(AM_CXXFLAGS) $(CXXFLAGS) \
	-O2 -c $< -fPIC -DPIC -o $@
endif

CLEANFILES = *_dict.cc *_dict.h *~ jetrec/*~ table/*~
