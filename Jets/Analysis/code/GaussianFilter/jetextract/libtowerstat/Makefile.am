# -*- mode: makefile; -*-

lib_LTLIBRARIES = libtowerstat.la

#####################################################################

libtowerstat_la_SOURCES = towerstat.cc towerstat_dict.cc
libtowerstat_la_CPPFLAGS = -D_LARGEFILE64_SOURCE @OFFLINE_CPPFLAGS@ \
	-I@srcdir@ -I@srcdir@/../libjetbase \
	-I@srcdir@/../libjetevent -I@srcdir@/../libhalf \
	-I@srcdir@/../libjetstat
libtowerstat_la_LDFLAGS = -L../libjetevent @ROOT_LDFLAGS@
libtowerstat_la_LIBADD = -ljetevent @ROOT_LIBS@ \
	-lfun4all -lrecal -lCNT -lCrkPID -lert -lspin -ltec -llvl2

noinst_HEADERS = towerstat/towerstat.h towerstat_linkdef.h
include_HEADERS =

%_dict.cc:	towerstat/%.h %_linkdef.h
	$(ROOTCINT) -f $@ -c $(libtowerstat_la_CPPFLAGS) $^

CLEANFILES = *_dict.cc *_dict.h *~ towerstat/*~
