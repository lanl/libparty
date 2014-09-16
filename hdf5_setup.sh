#!/bin/sh

# Use either ${HOME} or ${HOME}/software
HDF5_VERSION=hdf5-1.8.13
HDF5_INSTALL_BASE=${HOME}/${HDF5_VERSION}
HDF5_BUILD_BASE=${HOME}/HDF5_BUILD

if test ! -d "${HDF5_BUILD_BASE}"; then
   mkdir ${HDF5_BUILD_BASE}
fi
cd ${HDF5_BUILD_BASE}
pwd
if test -d "${HDF5_INSTALL_BASE}"; then
   echo "Removing old hdf5 install version"
   rm -rf ${HDF5_INSTALL_BASE}
fi
if test -d "${HDF5_VERSION}"; then
   echo "Removing old hdf5 build version"
   rm -rf ${HDF5_VERSION}
fi
echo "Untarring hdf5 -- have patience"
tar xzf ${HOME}/${HDF5_VERSION}.tar.gz
echo "Starting build of hdf5"
cd ${HDF5_VERSION}

if test ! -d "${HDF5_INSTALL_BASE}"; then
   mkdir ${HDF5_INSTALL_BASE}
fi
#CC=`which mpicc` ./configure --enable-parallel --prefix=${HDF5_INSTALL_BASE}
./configure --prefix=${HDF5_INSTALL_BASE}
make
make install


