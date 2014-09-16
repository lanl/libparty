#!/bin/csh
module load compilers/intel/14.0.2
module load mpi/openmpi-1.8.1-intel_14.0.2
set LD_LIBRARY_PATH=/home/kinsey/software/boost_1_54_0/stage/lib/:$LD_LIBRARY_PATH
source /projects/opt/intel/compilers/parallel_studio_xe_2013/composer_xe_2013_sp1.2.144/tbb/bin/tbbvars.csh intel64
