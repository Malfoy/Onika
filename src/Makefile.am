bin_PROGRAMS = onika

common_SOURCES = \
  common.h \
  zstr.hpp \
  strict_fstream.hpp \
  zstr_make_unique_polyfill.h \
  genome.h genome.cpp \
  onika_index.h onika_index.cpp

onika_SOURCES = \
  $(common_SOURCES) \
  optionparser.h \
  onika.cpp

MOSTLYCLEANFILES=*~
CLEANFILES=*~ $(noinst_SCRIPTS)
DISTCLEANFILES=*~
MAINTAINERCLEANFILES=*~

