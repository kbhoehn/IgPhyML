#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
# Original configurations provided by Dr. Stéphane Guindon. 
# Modifications done by Marcelo S Zanetti.

AC_PREREQ(2.61)
AC_INIT([IgPHYML],[1.0 201703.16],[kenneth.hoehn@gmail.com])
dnl AC_CONFIG_SRCDIR([simu.c])
AC_CONFIG_SRCDIR([src/simu.c])
AC_CONFIG_HEADER([config.h])

AM_INIT_AUTOMAKE
AC_CANONICAL_HOST

# Checks for programs.
AC_PROG_CC

# Checks for libraries.
AC_CHECK_LIB([m], [log])

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([float.h stdlib.h string.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_FUNC_SETVBUF_REVERSED
AC_FUNC_VPRINTF
AC_CHECK_FUNCS([floor pow rint sqrt strchr strstr])

# YAML
AC_ARG_WITH([yaml],[  --with-yaml    enable YAML])
AM_CONDITIONAL([USE_YAML], [test "$with_yaml" == yes]) 
if test "$with_yaml" == "yes"
then
    # Get library and include locations
    AC_ARG_WITH([yaml-path],
                  [AS_HELP_STRING([--with-yaml-path],
                                    [location of the YAML installation, defaults to /usr/local])],
                                    [YAML_FLAGS="$withval"],
                                    [YAML_FLAGS='/usr/local'])
    AC_SUBST([YAML_FLAGS]) 
    AC_DEFINE([USEYAML], [1], [Define if yaml is used.])
fi
# Added for cross-compilation (thanks to Alban Mancheron)
MACOS_BASE="/Developer"
AC_ARG_WITH([macos_base],
            [AS_HELP_STRING([--with-macos-base@<:@=MACOS_BASE@:>@],
                            Path to the MACOS base directory @<:@/Developer@:>@)],
            [MACOS_BASE="$withval"],
            [])
 
MACOS_SDK="$MACOS_BASE/SDKs/MacOSX10.4u.sdk/"
AC_ARG_WITH([macos_sdk],
            [AS_HELP_STRING([--with-macos-sdk@<:@=MACOS_SDK@:>@],
                            Path to the MACOS SDKs directory @<:@$(MACOS_BASE)/SDKs/MacOSX10.4u.sdk/@:>@)],
            [MACOS_SDK="$withval"],
            [])



case $host_os in
  darwin*)  CFLAGS+=" -I$MACOS_SDK/usr/include/gcc/                                    \
                      -I$MACOS_BASE/include/                                           \
                      -arch i386 -arch ppc                                             \
                      -Wno-long-long -Wno-variadic-macros"
            AC_DEFINE([MACOS], [1], [Running on a mac?])
            AM_CONDITIONAL([MACOS], [test 1 = 1]);;
dnl             LDFLAGS+=" -arch i386 -arch ppc";;
  mingw32*) CFLAGS+=" -std=gnu99"
            AM_CONDITIONAL([MACOS], [test 1 = 0]);;
  linux*|*) AM_CONDITIONAL([MACOS], [test 1 = 0]);;
esac

# End of cross-compilation macros

# CFLAGS 
# -g                   = debugging
# -pg                  = profiling
# -fopenmp             = openMP support 
# -O3                  = optimization level (fast option) 
# -Wfloat-equal        = Warn about = and <> comparisons of real numbers. 
# -fomit-frame-pointer = Don't keep the frame pointer in a register for functions that don't need one. This avoids the instructions to save, set up and restore frame pointers; it also makes an extra register available in many functions. It also makes debugging impossible on some machines.
# -funroll-loops       = Unroll loops whose number of iterations can be determined at compile time or upon entry to the loop. -funroll-loops implies -frerun-cse-after-loop. It also turns on complete loop peeling (i.e. complete removal of loops with small constant number of iterations). This option makes code larger, and may or may not make it run faster. 

# Recommended settings
CFLAGS="-O3 -fomit-frame-pointer -funroll-loops -fopenmp"
dnl CFLAGS="-O3 -fomit-frame-pointer -funroll-loops -static"
dnl CFLAGS = "-g"
dnl CFLAGS="-O3 -pg -funroll-loops"

AC_ARG_ENABLE([windows],[  --enable-windows  Compile windows],[windows=yes])
AM_CONDITIONAL([WANT_WINDOWS], [test "$windows" = yes]) 

AC_ARG_ENABLE([blas_omp],[  --enable-blas_omp    Compile blas_omp],[blas_omp=yes])
AM_CONDITIONAL([WANT_BLAS_OMP], [test "$blas_omp" = yes]) 
if test "$blas_omp" == "yes"
then
    AC_DEFINE([BLAS_OMP], [1], [Define if BLAS and OMP is used.])
fi

AC_ARG_ENABLE([omp],[  --enable-omp    Compile omp],[omp=yes])
AM_CONDITIONAL([WANT_OMP], [test "$omp" = yes]) 
if test "$omp" == "yes"
then
    AC_DEFINE([OMP], [1], [Define if OMP is used.])
fi

AC_ARG_ENABLE([blas],[  --enable-blas    Compile blas],[blas=yes])
AM_CONDITIONAL([WANT_BLAS], [test "$blas" = yes]) 
if test "$blas" == "yes"
then
    AC_DEFINE([BLAS], [1], [Define if BLAS is used.])
fi

AC_ARG_ENABLE([tiporder],[  --enable-tiporder    Compile tiporder],[tiporder=yes])
AM_CONDITIONAL([WANT_TIPORDER], [test "$tiporder" = yes]) 

AC_ARG_ENABLE([times],[  --enable-times    Compile PhyTime],[times=yes])
AM_CONDITIONAL([WANT_TIMES], [test "$times" = yes]) 

AC_ARG_ENABLE([part],[  --enable-part    Compile part],[part=yes])
AM_CONDITIONAL([WANT_PART], [test "$part" = yes]) 

AC_CONFIG_FILES([Makefile src/Makefile])
AC_OUTPUT
