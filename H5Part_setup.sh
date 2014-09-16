#!/bin/sh

# Get dist file from https://codeforge.lbl.gov/projects/h5part 

# Use either ${HOME} or ${HOME}/software
H5Part_VERSION=H5Part-1.6.6
H5Part_INSTALL_BASE=${HOME}/${H5Part_VERSION}
H5Part_BUILD_BASE=${HOME}/H5Part_BUILD

if test ! -d "${H5Part_BUILD_BASE}"; then
   mkdir ${H5Part_BUILD_BASE}
fi
cd ${H5Part_BUILD_BASE}
pwd
if test -d "${H5Part_INSTALL_BASE}"; then
   echo "Removing old H5Part install version"
   rm -rf ${H5Part_INSTALL_BASE}
fi
if test -d "${H5Part_VERSION}"; then
   echo "Removing old H5Part build version"
   rm -rf ${H5Part_VERSION}
fi
echo "Untarring H5Part -- have patience"
tar xzf ${HOME}/${H5Part_VERSION}.tar.gz
echo "Starting build of H5Part"
cd ${H5Part_VERSION}

if test ! -d "${H5Part_INSTALL_BASE}"; then
   mkdir ${H5Part_INSTALL_BASE}
fi
#CC=`which mpicc` CXX=`which mpicxx` MPICC=`which mpicc` MPICXX=`which mpicxx` ./configure --enable-parallel --prefix=${H5Part_INSTALL_BASE}
./configure  --enable-fortran --enable-shared --prefix=${H5Part_INSTALL_BASE}
make
make install


