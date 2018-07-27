# -*- mode: makefile; -*-

lib_LTLIBRARIES = libjetevent.la

#####################################################################

libjetevent_la_SOURCES = event.cc eventrec.cc event_dict.cc lorentz.cc \
	lorentz_dict.cc particle.cc particle_dict.cc snowmass.cc \
	snowmass_dict.cc #xefio.cc
libjetevent_la_CPPFLAGS = -I@srcdir@ -I@srcdir@/../libjetbase \
	-I@srcdir@/../libppmd -I@srcdir@/../libhalf
libjetevent_la_LDFLAGS = @ROOT_LDFLAGS@
libjetevent_la_LIBADD = @ROOT_LIBS@

noinst_HEADERS =
nobase_include_HEADERS = jetevent/event.h jetevent/lorentz.h \
	jetevent/particle.h jetevent/snowmass.h #jetevent/xefio.h
pkgdata_DATA =

%_dict.cc:	jetevent/%.h %_linkdef.h
	$(ROOTCINT) -f $@ -c $(libjetevent_la_CPPFLAGS) $^

CLEANFILES = *_dict.cc *_dict.h *~ jetevent/*~ table/*~
