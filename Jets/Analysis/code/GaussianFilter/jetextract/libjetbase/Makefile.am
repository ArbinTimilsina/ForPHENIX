# -*- mode: makefile; -*-

lib_LTLIBRARIES = libjetbase.la

#####################################################################

libjetbase_la_SOURCES = geometry_dict.cc hash.cc mem.cc num.cc \
	num_dict.cc specfunc.cc
libjetbase_la_CPPFLAGS = -I@srcdir@
libjetbase_la_LDFLAGS = @ROOT_LDFLAGS@
libjetbase_la_LIBADD = @ROOT_LIBS@

noinst_HEADERS =
nobase_include_HEADERS = jetbase/dbc.h jetbase/geometry.h jetbase/hash.h \
	jetbase/mem.h jetbase/num.h jetbase/specfunc.h
pkgdata_DATA =

%_dict.cc:	jetbase/%.h %_linkdef.h
	$(ROOTCINT) -f $@ -c $(libjetbase_la_CPPFLAGS) $^

CLEANFILES = *_dict.cc *_dict.h *~ jetbase/*~ table/*~
