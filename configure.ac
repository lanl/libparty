AC_INIT([libparty], [v1.0])

AC_CONFIG_AUX_DIR(./config)

dnl AUTOMAKE_OPTIONS([foreign])
AC_CONFIG_SRCDIR([cc_test.c])
dnl AM_INIT_AUTOMAKE([libparty], [1.0])
AM_INIT_AUTOMAKE([foreign])

AC_CONFIG_HEADERS([src/config.h config.h])

dnl
dnl Caching is usually WRONG for systems with cross-mounted file systems
dnl  the cache file may correspond to a different system.
dnl
define([AC_CACHE_LOAD], )dnl
define([AC_CACHE_SAVE], )dnl
dnl

#
echo "Configuring libparty with $ac_configure_args"
#

AC_ARG_ENABLE( debug,
 [  --enable-debug          compiles code with debug flags ],
 )

AC_ARG_ENABLE( bounds-check,
 [  --enable-bounds-check   compiles code with bounds check flags ],
 )

AC_ARG_ENABLE( prof,
 [  --enable-prof           compiles code with profiling flags ],
 )

AC_ARG_ENABLE( openmp,
 [  --enable-openmp         compiles code with openmp ],
 )

AC_ARG_ENABLE( mic-offload,
 [  --enable-mic-offload    compiles code with mic offload],
 )

AC_ARG_ENABLE( mic,
 [  --enable-mic    compiles code with mic ],
 )

AC_ARG_ENABLE( mpi, 
 [  --enable-mpi            compiles code with mpi library ],
 enable_mpi=no)

AC_ARG_ENABLE( opencl, 
 [  --enable-opencl         compiles code with opencl library ],
 enable_opencl=no)

AC_ARG_ENABLE( boost, 
 [  --enable-boost          compiles code with boost library ],
 )

dnl AC_ARG_ENABLE( hdf5, 
dnl  [  --enable-hdf5           compiles code with hdf5 library ],
dnl  enable_hdf5="yes"
dnl  )

AC_ARG_ENABLE( debug,
 [  --enable-scaling_tests  compiles code for scaling tests ],
 )

AC_ARG_ENABLE( H5Part,
 [  --enable-H5Part         compiles code with H5Part library ],
 enable_H5Part="yes"
 )

AC_ARG_ENABLE( detailed-timing, 
 [  --enable-detailed-timing compiles code with detailed timing on ],
 )

AC_ARG_ENABLE( write-mem-usage, 
 [  --enable-write-mem-usage compiles code with write memory usage on ],
 )

if test "${enable_mpi}" = "yes"; then
   enable_boost="yes"
fi

if test "${enable_mic_offload}" = "yes"; then
   enable_openmp="yes"
   enable_boost="no"
   AC_DEFINE(USE_MIC,1,"use mic hardware")
fi

if test "${enable_mic}" = "yes"; then
   enable_openmp="yes"
   enable_boost="no"
   AC_DEFINE(USE_MIC,1,"use mic hardware")
fi

if test "${enable_scaling_tests}" = "yes"; then
   AC_DEFINE(USE_SCALING_TESTS,1,"turns on scaling tests")
fi

#AM_WITH_DMALLOC

dnl so we don't have shared libraries
dnl AC_ENABLE_STATIC(yes)

# Checks for programs
AC_PROG_CC
AC_PROG_CXX
AC_PROG_FC
AC_C_INLINE

# Setting OpenMP compiler flags for C and C++
# requires configure 2.62
AC_LANG(C)
AC_OPENMP
AC_LANG_PUSH(C)
AC_LANG(C++)
AC_OPENMP
AC_LANG_POP(C++)

#AM_PROG_CC_C_O
AC_PROG_INSTALL
AC_PROG_MAKE_SET
AC_PROG_LIBTOOL
AC_PROG_RANLIB

AC_CONFIG_MACRO_DIR([m4])

INDENT="indent"
CPREPROCESSOR="cpp -P"

AC_CHECK_PROG([HAVE_ADDR2LINE], [addr2line], [yes])
if test "x$HAVE_ADDR2LINE" == "xyes"; then
   AC_DEFINE([HAVE_ADDR2LINE],[1],[Has addr2line program])
fi

# Checks for header files
AC_HEADER_STDC
AC_HEADER_STDBOOL

# Checks for typedefs, structures, and compiler characteristics
AC_SYS_LARGEFILE
AC_TYPE_SIZE_T

# Checks for library functions
AC_CHECK_FUNCS([floor])
AC_CHECK_FUNCS([pow])
AC_CHECK_FUNCS([sqrt])


AC_CANONICAL_HOST
case "$host" in
    *linux*)
dnl     if test "$CC" = "gcc"; then
dnl        if test -z "$F90"; then
dnl           F90="gfortran"
dnl           FCC="gfortran"
dnl           CXXFLAGS="-std=c++11 ${CXXFLAGS}"
dnl        fi
dnl     fi
dnl     if test "$CC" = "icc"; then
dnl        if test -z "$F90"; then
dnl           F90="ifort"
dnl        fi
dnl     fi
dnl     if test "$CC" = "mpicc"; then
dnl        if test -z "$F90"; then
dnl           F90="mpifort"
dnl        fi
dnl     fi
        dnl echo "CXX ",$CXX
        dnl echo "CC ",$CC
        dnl echo "FCC ",$FCC
        dnl echo "FCFLAGS ",$FCFLAGS
        dnl echo "F90 ",$F90
        if test "${OPENMP_CFLAGS}x" != "x"; then
           CFLAGS="${OPENMP_CFLAGS} ${CFLAGS}"
           CXXFLAGS="${OPENMP_CXXFLAGS} ${CXXFLAGS}"
        fi
           
        if test "`basename ${CC}`" = "icc"; then
           CFLAGS="-wd1478 ${CFLAGS}"
           LIBS="-L/projects/opt/intel/compilers/parallel_studio_xe_2013/composer_xe_2013_sp1.2.144/tbb/lib/intel64/gcc4.4 ${LIBS}"
           if test "${RPATH_LIBS}x" = "x"; then
             RPATH_LIBS="/projects/opt/intel/compilers/parallel_studio_xe_2013/composer_xe_2013_sp1.2.144/tbb/lib/intel64/gcc4.4"
           else
             RPATH_LIBS="/projects/opt/intel/compilers/parallel_studio_xe_2013/composer_xe_2013_sp1.2.144/tbb/lib/intel64/gcc4.4:${RPATH_LIBS}"
           fi
           TBBROOT="/projects/opt/intel/compilers/parallel_studio_xe_2013/composer_xe_2013_sp1.2.144/tbb"
           MIC_LD_LIBRARY_PATH="/projects/opt/intel/compilers/parallel_studio_xe_2013/composer_xe_2013_sp1.2.144/tbb/lib/mic"
        fi
        if test "`basename ${CXX}`" = "icpc"; then
           CXXFLAGS="-wd1478 ${CXXFLAGS}"
        fi

        CFLAGS="-static ${CFLAGS} -tbb"
        CXXFLAGS="-static -std=c++0x -ansi -tbb ${CXXFLAGS}"
        dnl CPPFLAGS="-I${INCDIR}"
        dnl FCFLAGS="-static -fPIC ${FCFLAGS} -tbb"

        if test "${enable_debug}" = "yes"; then
          CFLAGS="-g -DDEBUG=1 -Wall ${CFLAGS}"
          CXXFLAGS="-g -DDEBUG=1 -Wall ${CXXFLAGS}"
          if test "$CC" = "gcc"; then
             FCFLAGS="-g ${FCFLAGS}"
          fi
        else
          CFLAGS="-g -O2 ${CFLAGS}"
          CXXFLAGS="-g -O2 ${CXXFLAGS}"
          if test "$CC" = "gcc"; then
             FCFLAGS="-g -O2 ${FCFLAGS}"
          else
             FCFLAGS="-O2 ${FCFLAGS}"
          fi
        fi
        if test "${enable_prof}" = "yes"; then
          CFLAGS="-pg ${CFLAGS}"
          CXXFLAGS="-pg ${CXXFLAGS}"
          FCFLAGS="-pg ${FCFLAGS}"
          AM_LDFLAGS="-pg ${LFLAGS}"
        fi
        if test "${enable_openmp}" = "yes"; then
          CFLAGS="-fopenmp ${CFLAGS}"
          CXXFLAGS="-fopenmp ${CXXFLAGS}"
          FCFLAGS="-openmp ${FCFLAGS}"
          LDFLAGS="-openmp ${LDFLAGS}"
        fi

        if test "${enable_bounds_check}" = "yes"; then
          CFLAGS="-fbounds-check ${CFLAGS}"
          CXXFLAGS="-fbounds-check ${CXXFLAGS}"
        fi

        RPATH_PREFIX="-Wl,"

        ;;
    *)
        echo "Unknown machine"
        exit 1
        ;;
esac

if test "${enable_mic_offload}" = "yes"; then
   CFLAGS="-offload-attribute-target=mic ${CFLAGS}"
   dnl CXXFLAGS="-offload-attribute-target=mic ${CXXFLAGS}"
   dnl FCFLAGS="-offload-attribute-target=mic ${FCFLAGS}"
fi

if test "${enable_mic}" = "yes"; then
   CFLAGS="-mmic ${CFLAGS}"
   dnl CXXFLAGS="-mmic ${CXXFLAGS}"
   dnl FCFLAGS="-mmic ${FCFLAGS}"
fi

if test "${enable_detailed_timing}" = "yes"; then
   CFLAGS="-DDETAILED_TIMING=1 ${CFLAGS}"
fi

if test "${enable_write_mem_usage}" = "no"; then
   CFLAGS="-DWRITE_MEM_USAGE=0 ${CFLAGS}"
fi

#if test "${enable_opencl}" != "no"; then
if test "${enable_opencl}" == "yes"; then
  case "$host" in
    *linux*)
      found="no"
      AC_CHECK_PROG(found, lspci, yes)
      if test ! -z "${OPENCL_INCLUDE}" ; then 
        AM_CPPFLAGS="${AM_CPPFLAGS} -I${OPENCL_INCLUDE}"
        if test ! -z "${OPENCL_LIBS}" ; then
          LIBS="${LIBS} -L${OPENCL_LIBS} -lOpenCL"
          if test "${RPATH_LIBS}x" = "x"; then
            RPATH_LIBS="${OPENCL_LIBS}"
          else
            RPATH_LIBS="${RPATH_LIBS}:${OPENCL_LIBS}"
          fi
        fi
      elif test "$found" = "no"; then
         found="no";
         AC_CHECK_FILE(/projects/opt/cuda/cuda-5.5/include/CL/cl.h, AM_CPPFLAGS="${AM_CPPFLAGS} -I/projects/opt/cuda/cuda-5.5/include" && found="yes")
         if test "$found" = "no"; then
            AC_CHECK_FILE( /projects/opt/opencl1d2_additional_files/include/CL/cl.h, AM_CPPFLAGS="${AM_CPPFLAGS} -I/projects/opt/opencl1d2_additional_files/include" && found="yes")
         fi
         if test "$found" = "no"; then
            AC_CHECK_FILE( /projects/opt/AMD_APP_SDK/AMDAPP-2.8/include/CL/cl.h, AM_CPPFLAGS="${AM_CPPFLAGS} -I/projects/opt/AMD_APP_SDK/AMDAPP-2.8/include" && found="yes")
         fi
         if test "$found" = "no"; then
            AC_CHECK_FILE( /opt/intel/opencl-1.2-3.0.67279/include/CL/cl.h, AM_CPPFLAGS="${AM_CPPFLAGS} -I/opt/intel/opencl-1.2-3.0.67279/include" && found="yes")
         fi
         found="no"
         if test "$found" = "no"; then
           AC_CHECK_FILE(/usr/lib/libOpenCL.so, LIBS="${LIBS} -lOpenCL" && found="yes")
           if test "${found}" = "yes"; then
             if test "${RPATH_LIBS}x" = "x"; then
               RPATH_LIBS="/usr/lib"
             else
               RPATH_LIBS="${RPATH_LIBS}:/usr/lib"
             fi
           fi
         fi
         if test "$found" = "no"; then
           AC_CHECK_FILE(/usr/lib64/libOpenCL.so, LIBS="${LIBS} -lOpenCL" && found="yes")
           if test "${found}" = "yes"; then
             if test "${RPATH_LIBS}x" = "x"; then
               RPATH_LIBS="/usr/lib64"
             else
               RPATH_LIBS="${RPATH_LIBS}:/usr/lib64"
             fi
           fi
         fi
         if test "$found" = "no"; then
           AC_CHECK_FILE(/projects/opt/AMD_APP_SDK/AMDAPP-2.8/lib/x86_64/libOpenCL.so, LIBS="${LIBS} -L/projects/opt/AMD_APP_SDK/AMDAPP-2.8/lib/x86_64 -lOpenCL" && found="yes")
           if test "${found}" = "yes"; then
             if test "${RPATH_LIBS}x" = "x"; then
               RPATH_LIBS="/projects/opt/AMD_APP_SDK/AMDAPP-2.8/lib/x86_64"
             else
               RPATH_LIBS="${RPATH_LIBS}:/projects/opt/AMD_APP_SDK/AMDAPP-2.8/lib/x86_64"
             fi
           fi
         fi
         if test "$found" = "no"; then
           AC_CHECK_FILE(/opt/intel/opencl-1.2-3.0.67279/lib64/libOpenCL.so, LIBS="${LIBS} -L/opt/intel/opencl-1.2-3.0.67279/lib64 -L/opt/intel/opencl-1.2-3.0.67279/libmic" && found="yes")
           if test "${found}" = "yes"; then
             if test "${RPATH_LIBS}x" = "x"; then
               RPATH_LIBS="/opt/intel/opencl-1.2-3.0.67279/lib64"
             else
               RPATH_LIBS="${RPATH_LIBS}:/opt/intel/opencl-1.2-3.0.67279/lib64"
             fi
           fi
         fi
      elif test "`lspci | grep "VGA compatible controller" | grep -i nVidia | wc -l`" != "0"; then
         found="no"
         AC_CHECK_FILE(/usr/local/NVIDIA_GPU_Computing_SDK/OpenCL/common/inc/CL/cl.h, AM_CPPFLAGS="${AM_CPPFLAGS} -I/usr/local/NVIDIA_GPU_Computing_SDK/OpenCL/common/inc" && found="yes")
         if test "$found" = "no"; then
           AC_CHECK_FILE(/usr/local/cuda/include/CL/cl.h, AM_CPPFLAGS="${AM_CPPFLAGS} -I/usr/local/cuda/include" && found="yes")
         fi
         if test "$found" = "no"; then
           AC_CHECK_FILE(/projects/opt/cuda/cuda-5.5/include/CL/cl.h, AM_CPPFLAGS="${AM_CPPFLAGS} -I/projects/opt/cuda/cuda-5.5/include
" && found="yes")
	fi
	if test "$found" = "no"; then
           AC_CHECK_FILE( /projects/opt/opencl1d2_additional_files/include/CL/cl.h, AM_CPPFLAGS="${AM_CPPFLAGS} -I/projects/opt/opencl1d2_additional_files/include" && found="yes")
         fi
         if test "$found" = "no"; then
           AC_CHECK_FILE(/home/cuda/sdk4.0.13/OpenCL/common/inc/CL/cl.h, AM_CPPFLAGS="${AM_CPPFLAGS} -I/home/cuda/sdk4.0.13/OpenCL/common/inc" && found="yes")
         fi
         found="no"
         if test "$found" = "no"; then
           AC_CHECK_FILE(/usr/lib/libOpenCL.so, LIBS="${LIBS} -lOpenCL" && found="yes")
           if test "${found}" = "yes"; then
             if test "${RPATH_LIBS}x" = "x"; then
               RPATH_LIBS="/usr/lib"
             else
               RPATH_LIBS="${RPATH_LIBS}:/usr/lib"
             fi
           fi
         fi
         if test "$found" = "no"; then
           AC_CHECK_FILE(/usr/lib64/libOpenCL.so, LIBS="${LIBS} -lOpenCL" && found="yes")
           if test "${found}" = "yes"; then
             if test "${RPATH_LIBS}x" = "x"; then
               RPATH_LIBS="/usr/lib64"
             else
               RPATH_LIBS="${RPATH_LIBS}:/usr/lib64"
             fi
           fi
         fi
      elif test "`lspci | grep "VGA compatible controller" | grep ATI | wc -l`" != "0"; then 
         found="no"
         if test "$found" = "no"; then
           AC_CHECK_FILE(/home/opt/AMDAPP_v2.6/include/CL/cl.h, AM_CPPFLAGS="${AM_CPPFLAGS} -I/home/opt/AMDAPP_v2.6/include" && found="yes")
         fi
         if test "$found" = "no"; then
           AC_CHECK_FILE(/projects/opt/AMD_APP_SDK/AMDAPP_v2.6/include/CL/cl.h, AM_CPPFLAGS="${AM_CPPFLAGS} -I/projects/opt/AMD_APP_SDK/AMDAPP_v2.6/include" && found="yes")
         fi
         found="no"
         if test "$found" = "no"; then
           AC_CHECK_FILE(/home/opt/AMDAPP_v2.6/lib/x86_64/libOpenCL.so, LIBS="${LIBS} -L/home/opt/AMDAPP_v2.6/lib/x86_64 -lOpenCL" && found="yes")
           if test "${found}" = "yes"; then
             if test "${RPATH_LIBS}x" = "x"; then
               RPATH_LIBS="/home/opt/AMDAPP_v2.6/lib/x86_64"
             else
               RPATH_LIBS="${RPATH_LIBS}:/home/opt/AMDAPP_v2.6/lib/x86_64"
             fi
           fi
         fi
         if test "$found" = "no"; then
           AC_CHECK_FILE(/projects/opt/AMD_APP_SDK/AMDAPP_v2.6/lib/x86_64/libOpenCL.so, LIBS="${LIBS} -L/projects/opt/AMD_APP_SDK/AMDAPP_v2.6/lib/x86_64 -lOpenCL" && found="yes")
           if test "${found}" = "yes"; then
             if test "${RPATH_LIBS}x" = "x"; then
               RPATH_LIBS="/projects/opt/AMD_APP_SDK/AMDAPP_v2.6/lib/x86_64"
             else
               RPATH_LIBS="${RPATH_LIBS}:/projects/opt/AMD_APP_SDK/AMDAPP_v2.6/lib/x86_64"
             fi
           fi
         fi
      fi
      ;;
    *darwin*)
      LIBS="${LIBS} -framework OpenCL"
      ;;
  esac


dnl ---------------------------------------------------------------
   AC_MSG_CHECKING(for OpenCL double type)

   SAVE_C_FLAGS=$CFLAGS
   CFLAGS="$CFLAGS $AM_CPPFLAGS"

   AC_LANG([C])

   AC_COMPILE_IFELSE([AC_LANG_PROGRAM(
      [
       [#ifdef __APPLE_CC__]
       [#include <OpenCL/OpenCL.h>]
       [#else]
       [#include <CL/cl.h>]
       [#endif]
      ],
      [[cl_double x;]])],
      ac_cv_opencl_double=yes, ac_cv_opencl_double=no
   )
   if test "x$ac_cv_opencl_double" = "xyes"; then
      AC_DEFINE(HAVE_CL_DOUBLE,1,[Has OpenCL double type])
   else
      AC_DEFINE(NO_CL_DOUBLE,1,[Does not have OpenCL double type])
   fi

   CFLAGS=$SAVE_C_FLAGS

   AC_MSG_RESULT($ac_cv_opencl_double)
dnl ---------------------------------------------------------------
fi

if test "${enable_mpi}" != "no"; then
dnl ---------------------------------------------------------------
   AC_PATH_PROG(MPI_BIN, mpicc)

  if test ! -z "${MPI_BIN}"; then
    MPIHOME=`dirname ${MPI_BIN}`
    if test "${host}" = "darwin"; then
       MPIHOME=`dirname ${MPIHOME}`
    fi
  fi

  dnl Find mpi include
  dnl ----------------
  FOUND_MPI_INCLUDE="no"
  AC_CHECK_FILE(${MPIHOME}/include/mpi.h, INCLUDES="${INCLUDES} -I${MPIHOME}/include" && FOUND_MPI_INCLUDE="yes")
  if test "${FOUND_MPI_INCLUDE}" = "no"; then
    AC_CHECK_FILE(${MPI_ROOT}/include/mpi.h, MPIHOME=${MPI_ROOT} && INCLUDES="${INCLUDES} -I${MPIHOME}/include" && FOUND_MPI_INCLUDE="yes")
  fi

  dnl Find mpi c++ library
  dnl ----------------
  FOUND_MPI_CPP_LIB="no"
  AC_CHECK_FILE(${MPIHOME}/lib/libmpi_cxx.so, MPI_LIBS="${MPI_LIBS} -L${MPIHOME}/lib -lmpi_cxx" && FOUND_MPI_CPP_LIB="yes" && MPI_LIB_LOCATE="${MPIHOME}")
  if test "${FOUND_MPI_CPP_LIB}" = "no"; then
     AC_CHECK_FILE(${MPIHOME}/lib/libmpi_cxx.dylib, MPI_LIBS="${MPI_LIBS} -L${MPIHOME}/lib -lmpi_cxx" && FOUND_MPI__CPP_LIB="yes" && MPI_LIB_LOCATE="${MPIHOME}")
  fi

  dnl Find mpi library
  dnl ----------------
  FOUND_MPI_LIB="no"
  if test "${MPI_LIB_LOCATE}x" = "${MPIHOME}x"; then
     AC_CHECK_FILE(${MPIHOME}/lib/libmpi.so, MPI_LIBS="${MPI_LIBS} -lmpi" && FOUND_MPI_LIB="yes")
     if test "${FOUND_MPI_LIB}" = "no"; then
        AC_CHECK_FILE(${MPIHOME}/lib/libmpi.dylib, MPI_LIBS="${MPI_LIBS} -lmpi" && FOUND_MPI_LIB="yes")
     fi
  else
     AC_CHECK_FILE(${MPIHOME}/lib/libmpi.so, MPI_LIBS="${MPI_LIBS} -L${MPIHOME}/lib -lmpi" && FOUND_MPI_LIB="yes")
     if test "${FOUND_MPI_LIB}" = "no"; then
        AC_CHECK_FILE(${MPIHOME}/lib/libmpi.dylib, MPI_LIBS="${MPI_LIBS} -L${MPIHOME}/lib -lmpi" && FOUND_MPI_LIB="yes")
     fi
  fi

 
   AC_MSG_CHECKING(MPI compile)

   SAVE_C_FLAGS=$CFLAGS
   CFLAGS="$CFLAGS $INCLUDES"

   AC_LANG([C])

   AC_COMPILE_IFELSE([AC_LANG_PROGRAM(
      [
       [#include <mpi.h>]
      ],
      [
       [int argc;]
       [char **argv;]
       [MPI_Init(&argc, &argv);]
      ])], 
      ac_cv_mpi=yes, ac_cv_mpi=no
   )
   if test "x$ac_cv_mpi" = "xyes"; then
      dnl AC_DEFINE(HAVE_MPI,1,[Has MPI libraries])
      if test "${FOUND_MPI_LIB}" = "no"; then
         MPI_LIBS="${MPI_LIBS} -lmpi_cxx -lmpi"
      fi
   fi

   if test "${ac_cv_mpi}" == "no"; then
      CFLAGS=$SAVE_C_FLAGS
   fi

   if test "${ac_cv_mpi}" == "yes"; then
      CXXFLAGS="$CXXFLAGS $INCLUDES"
   fi

   AC_MSG_RESULT($ac_cv_mpi)
dnl ---------------------------------------------------------------
fi

AC_MSG_NOTICE("enable boost is ${enable_boost})
if test "${enable_boost}" != "no"; then
   FOUND_BOOST_LIB="no"
   BOOST_LIB_PATH_BASE="${HOME}/boost_1_54_0"
   AC_CHECK_FILE(${BOOST_LIB_PATH_BASE}/stage/lib/libboost_mpi.a,
      BOOST_LIB_PATH="${BOOST_LIB_PATH_BASE}/stage/lib"
      BOOST_INC_PATH="${BOOST_LIB_PATH_BASE}"
      FOUND_BOOST_LIB="yes")
   if test "${FOUND_BOOST_LIB}" = "no"; then
      BOOST_LIB_PATH_BASE="${HOME}/software/boost_1_54_0"
      AC_CHECK_FILE(${BOOST_LIB_PATH_BASE}/stage/lib/libboost_mpi.a,
         BOOST_LIB_PATH="${BOOST_LIB_PATH_BASE}/stage/lib"
         BOOST_INC_PATH="${BOOST_LIB_PATH_BASE}"
         FOUND_BOOST_LIB="yes")
   fi
   if test "${FOUND_BOOST_LIB}" = "no"; then
      BOOST_LIB_PATH_BASE="/home/kinsey/software/boost_1_54_0"
      AC_CHECK_FILE(${BOOST_LIB_PATH_BASE}/stage/lib/libboost_mpi.a,
         BOOST_LIB_PATH="${BOOST_LIB_PATH_BASE}/stage/lib"
         BOOST_INC_PATH="${BOOST_LIB_PATH_BASE}"
         FOUND_BOOST_LIB="yes")
   fi

   if test ! -z "${BOOST_LIB_PATH}"; then
      LIBS="${LIBS} -L${BOOST_LIB_PATH} -lboost_mpi -lboost_serialization"
      CFLAGS="${CFLAGS} -I${BOOST_INC_PATH}"
      CXXFLAGS="${CXXFLAGS} -I${BOOST_INC_PATH}"
      if test "${RPATH_LIBS}x" = "x"; then
         RPATH_LIBS="${BOOST_LIB_PATH}"
      else
         RPATH_LIBS="${RPATH_LIBS}:${BOOST_LIB_PATH}"
      fi
   fi
fi

if test "${enable_hdf5}" != "no"; then
   FOUND_HDF5_LIB="no"
   if test "${FOUND_HDF5_LIB}" = "no"; then
      HDF5_LIB_PATH_BASE="${HOME}/hdf5-1.8.13"
      AC_CHECK_FILE(${HDF5_LIB_PATH_BASE}/lib/libhdf5.a,
         HDF5_LIB_PATH="${HDF5_LIB_PATH_BASE}/lib"
         HDF5_INC_PATH="${HDF5_LIB_PATH_BASE}/include"
         FOUND_HDF5_LIB="yes")
   fi
   if test "${FOUND_HDF5_LIB}" = "no"; then
      HDF5_LIB_PATH_BASE="${HOME}/software/hdf5-1.8.13"
      AC_CHECK_FILE(${HDF5_LIB_PATH_BASE}/lib/libhdf5.a,
         HDF5_LIB_PATH="${HDF5_LIB_PATH_BASE}/lib"
         HDF5_INC_PATH="${HDF5_LIB_PATH_BASE}/include"
         FOUND_HDF5_LIB="yes")
   fi
   if test "${FOUND_HDF5_LIB}" = "no"; then
      HDF5_LIB_PATH_BASE="/usr/lib64"
      AC_CHECK_FILE(/usr/lib64/libhdf5.so,
         dnl Do not need to add to compiler options -- in default locations
         dnl HDF5_LIB_PATH="/usr/lib64" 
         dnl HDF5_INC_PATH="/usr/include"
         FOUND_HDF5_LIB="yes")
   fi

   if test ! -z "${HDF5_LIB_PATH}"; then
      LIBS="${LIBS} -L${HDF5_LIB_PATH} -lhdf5"
      CFLAGS="${CFLAGS} -I${HDF5_INC_PATH}"
      CXXFLAGS="${CXXFLAGS} -I${HDF5_INC_PATH}"
   fi
fi

if test "${enable_H5Part}" != "no"; then
   FOUND_H5PART_LIB="no"
   H5PART_LIB_PATH_BASE="${HOME}/H5Part-1.6.6"
   AC_CHECK_FILE(${H5PART_LIB_PATH_BASE}/lib/libH5Part.a,
      H5PART_LIB_PATH="${H5PART_LIB_PATH_BASE}/lib"
      H5PART_INC_PATH="${H5PART_LIB_PATH_BASE}/include"
      FOUND_H5PART_LIB="yes")
   if test "${FOUND_H5PART_LIB}" = "no"; then
      H5PART_LIB_PATH_BASE="${HOME}/software/H5Part-1.6.6"
      AC_CHECK_FILE(${H5PART_LIB_PATH_BASE}/lib/libH5Part.a,
         H5PART_LIB_PATH="${H5PART_LIB_PATH_BASE}/lib"
         H5PART_INC_PATH="${H5PART_LIB_PATH_BASE}/include"
         FOUND_H5PART_LIB="yes")
   fi
   if test "${FOUND_H5PART_LIB}" = "no"; then
      H5PART_LIB_PATH_BASE="/home/kinsey/software/H5Part-1.6.6/build"
      AC_CHECK_FILE(${H5PART_LIB_PATH_BASE}/lib/libH5Part.a,
         H5PART_LIB_PATH="${H5PART_LIB_PATH_BASE}/lib"
         H5PART_INC_PATH="${H5PART_LIB_PATH_BASE}/include"
         FOUND_H5PART_LIB="yes")
   fi

   if test ! -z "${H5PART_LIB_PATH}"; then
      LIBS="${LIBS} -L${H5PART_LIB_PATH} -lH5Part"
      CFLAGS="${CFLAGS} -I${H5PART_INC_PATH}"
      CXXFLAGS="${CXXFLAGS} -I${H5PART_INC_PATH}"
      if test "${RPATH_LIBS}x" = "x"; then
         RPATH_LIBS="${H5PART_LIB_PATH}"
      else
         RPATH_LIBS="${RPATH_LIBS}:${H5PART_LIB_PATH}"
      fi
   else
      enable_H5Part="no"
   fi
   if test "${enable_H5Part}" != "no"; then
      AC_DEFINE(HAVE_H5PART,1,"Have H5Part library")
   fi
fi

if test "${RPATH_LIBS}x" = "x"; then
   RPATH_LIBS="`pwd`/lib"
else
   RPATH_LIBS="${RPATH_LIBS}:`pwd`/lib"
fi

if test ! -z "${RPATH_LIBS}"; then
   AC_MSG_NOTICE("RPATH IS ${RPATH_LIBS}")
   LIBS="${LIBS} ${RPATH_PREFIX}-rpath ${RPATH_PREFIX}${RPATH_LIBS}"
fi

LIBS="${LIBS} -lm -lz"

AC_SUBST(CC)
AC_SUBST(CPP)
AC_SUBST(F90)
AC_SUBST(CFLAGS)
dnl AC_SUBST(AM_CPPFLAGS)
AC_SUBST(CPPFLAGS)
AC_SUBST(F90FLAGS)
AC_SUBST(FCFLAGS)
AC_SUBST(LDFLAGS)
AC_SUBST(LIBS)
AC_SUBST(INDENT)
AC_SUBST(CPREPROCESSOR)

AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([src/Makefile])

AC_OUTPUT

