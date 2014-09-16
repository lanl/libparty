all: libparty cctest cpptest f90test mpitest

bin_PROGRAMS = cctest cpptest f90test mpitest

SUBDIRS = src

cctest_SOURCES = cc_test.c
cctest_LDADD = -Lsrc -lparty -L/usr/lib64

cpptest_SOURCES = cpp_test.cc
cpptest_LDADD = -Lsrc -lparty -L/usr/lib64

f90test_SOURCES = fort_test.f90
f90test_FCFLAGS = -Isrc
f90test_LDADD = -Lsrc -lparty -lstdc++ -ltbb
f90test_LDFLAGS = -openmp

mpitest_CXX = mpic++
mpitest_SOURCES = mpi_test.cc
mpitest_LDADD = -Lsrc -lparty -L/usr/lib64

libparty:
	cd src; make

H5Part:
	H5Part_setup.sh

boost:
	boost_setup.sh

clean-local:
	rm -f particles_soa.h5part
	rm -f cctest.out
	rm -f cpptest.out
	rm -f f90test.out
	rm -f mpitest.out

distclean-local:
	rm -rf lib include src/libparty.modmic
	rm -rf autom4te.cache config.h.in~
	rm -rf autom4te.cache

check:
	@rm -f particles_soa.h5part
	./cctest > cctest.out
	@if test "`diff cctest.out test/gold_cctest.out | wc -l`" = 0; then echo "Passed cctest"; else echo "Failed cctest";  fi
	./cpptest > cpptest.out
	@if test "`diff cpptest.out test/gold_cpptest.out | wc -l`" = 0; then echo "Passed cpptest"; else echo "Failed cpptest";  fi
	@rm -f particles_soa.h5part
	./f90test > f90test.out
	@if test "`diff f90test.out test/gold_f90test.out | wc -l`" = 0; then echo "Passed f90test"; else echo "Failed f90test";  fi
	@rm -f particles_soa.h5part
	mpirun -n 2 mpitest > mpitest.out
	@grep -v "^Reported" mpitest.out > mpitest.tmpout
	@if test "`diff mpitest.tmpout test/gold_mpitest.out | wc -l`" = 0; then echo "Passed mpitest"; else echo "Failed mpitest";  fi
	@rm -f particles_soa.h5part
