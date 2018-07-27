# -*- mode: makefile; -*-

lib_LTLIBRARIES = libextract.la

#####################################################################

libextract_la_SOURCES = extract.cc extract_dict.cc extract_first.cc \
	extract_first_dict.cc
libextract_la_CPPFLAGS = -D_LARGEFILE64_SOURCE @OFFLINE_CPPFLAGS@ \
	-I@srcdir@ -I@srcdir@/../libjetbase \
	-I@srcdir@/../libjetevent -I@srcdir@/../libjetrec \
	-I@srcdir@/../libjetstat
libextract_la_LDFLAGS = -L../libjetevent @ROOT_LDFLAGS@
libextract_la_LIBADD = -ljetevent @ROOT_LIBS@ \
	-lfun4all -lrecal -lCNT -lCrkPID -lert -lspin -ltec -llvl2

noinst_HEADERS = extract/extract.h extract/extract_first.h \
	extract/run_5_p_p_tower_map.h \
	extract/run_5_p_p_tower_scale.h extract_linkdef.h
include_HEADERS =

%_dict.cc:	extract/%.h %_linkdef.h
	$(ROOTCINT) -f $@ -c $(libextract_la_CPPFLAGS) $^

CLEANFILES = *_dict.cc *_dict.h *~ extract/*~
