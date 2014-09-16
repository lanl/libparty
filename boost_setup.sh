#!/bin/sh

# Use either ${HOME} or ${HOME}/software
BOOST_INSTALL_BASE=${HOME}
BOOST_VERSION=boost_1_54_0

if test ! -d "${BOOST_INSTALL_BASE}"; then
   mkdir ${BOOST_INSTALL_BASE}
fi
cd ${BOOST_INSTALL_BASE}
pwd
if test -d "${BOOST_VERSION}"; then
   echo "Removing old boost version"
fi
echo "Untarring boost -- have patience"
tar xzf ${BOOST_VERSION}.tar.gz
echo "Starting build of boost"
cd ${BOOST_VERSION}
echo "" >> tools/build/v2/user-config.jam
echo "using mpi ;" >> tools/build/v2/user-config.jam
CC=icc f90=ifort ./bootstrap.sh --with-toolset=intel-linux
CC=icc f90=ifort ./b2 toolset=intel --with-mpi
