ACLOCAL_AMFLAGS = -I m4
# -*- mode: makefile; -*-

INCLUDES = -I$(includedir) -I$(OFFLINE_MAIN)/include -I$(ROOTSYS)/include

lib_LTLIBRARIES = libJetAnalyzer.la

libJetAnalyzer_la_SOURCES = JetAnalyzer.C JetAnalyzer_dict.C

libJetAnalyzer_la_CPPFLAGS = -D_LARGEFILE64_SOURCE @OFFLINE_CPPFLAGS@

libJetAnalyzer_la_LDFLAGS = @ROOT_LDFLAGS@ @OFFLINE_LDFLAGS@

libJetAnalyzer_la_LIBADD = -ljetevent @ROOT_LIBS@

noinst_HEADERS = JetAnalyzer.h JetAnalyzerLinkDef.h JetRecoTrackQuality.h

include_HEADERS =

%_dict.C:	%.h %LinkDef.h
	$(ROOTCINT) -f $@ -c $(libJetAnalyzer_la_CPPFLAGS) $^

CLEANFILES = *_dict.C *_dict.h *~ *~

clean-master:
	rm -r *_dict.* Makefile.in config.* configure lib* ltmain.sh Makefile aclocal.m4 m4 depcomp *~ autom4te.cache  install-sh missing .deps .libs *root macro/*root macro/*~ stamp-h1