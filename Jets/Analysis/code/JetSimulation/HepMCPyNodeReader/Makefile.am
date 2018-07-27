AUTOMAKE_OPTIONS = foreign

AM_CXXFLAGS = -DG4INTY_USE_XT -DG4VIS_USE_OPENGL -DG4UI_USE_TCSH -DG4UI_USE_XM -DG4VIS_USE_OPENGLXM -DG4VIS_USE_RAYTRACERX -DG4VIS_USE_OPENGLX 

# List of shared libraries to produce
lib_LTLIBRARIES = \
    libHepMCPyNodeReader.la

INCLUDES = \
    -I$(includedir) \
    -I$(OFFLINE_MAIN)/include \
    -I$(ROOTSYS)/include

libHepMCPyNodeReader_la_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib

libHepMCPyNodeReader_la_LIBADD = \
  -lphool \
  -lSubsysReco \
  -lfun4all \
  -lPHPythia \
  -lPHHijing \
  -lphhepmc

libHepMCPyNodeReader_la_SOURCES = \
    HepMCPyNodeReader_Dict.C \
    HepMCPyNodeReader.C

##############################################
# please add new classes in alphabetical order

include_HEADERS = \
  HepMCPyNodeReader.h

################################################
# linking tests

noinst_PROGRAMS = \
  testexternals_HepMCPyNodeReader

testexternals_HepMCPyNodeReader_SOURCES = testexternals.C
testexternals_HepMCPyNodeReader_LDADD = libHepMCPyNodeReader.la

testexternals.C:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

##############################################
# please add new classes in alphabetical order

clean-local:
	rm -f *Dict* $(BUILT_SOURCES)

# rules to make root dictionaries
%_Dict.C: %.h %LinkDef.h
	rootcint -f $@ -c $(DEFAULT_INCLUDES) $(INCLUDES) $^
