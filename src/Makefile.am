default: libparty.la
all: libparty.la

AM_MAKEFLAGS = -j 4

lib_LTLIBRARIES = libparty.la
include_HEADERS = libparty.h

F90_SRCS = libparty_bindings_soa.f90
CXX_SRCS = libparty.cc hashtable.cc libparty_wrappers.cc libparty_soa.cc libparty_h5part_soa.cc hashtable_soa.cc calcdens_soa.cc
H_SRCS = libparty.h libparty.hh hashutil.hh libparty_wrappers.h MersenneTwister.hh

libparty_la_SOURCES = ${F90_SRCS} ${CXX_SRCS} ${H_SRCS}

clean-local:
	rm -rf *.mod

DISTCLEANFILES = config.h.in~
