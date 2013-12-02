dnl -*-Autoconf-*-
dnl configure.ac for paradias
dnl Time-stamp: <2013-03-31 23:47:12 takeshi>
dnl Author: Takeshi NISHIMATSU

# Initialize.
AC_PREREQ(2.61)
AC_INIT([PARADIAS by Takeshi NISHIMATSU], 0.00.04, t-nissie{at}imr.tohoku.ac.jp, paradias)
AC_CONFIG_SRCDIR(src/paradias.F)
AM_INIT_AUTOMAKE
AM_CONFIG_HEADER(config.h)

# Some OS does not see /usr/local/lib.
LDFLAGS="$LDFLAGS -L/usr/local/lib"

# fftw3.f must be in /usr/include/, /usr/local/include/.
# Or specify the directory explicitly as:
# $ CPPFLAGS=-I/some/wherer ./configure
# AC_CHECK_HEADER(fftw3.f) does not work in Fortran.
CPPFLAGS="$CPPFLAGS -I/usr/include -I/usr/local/include"
AC_ARG_VAR(CPPFLAGS, [C preprocessor (CPP) Flags, e.g. -l<include dir>])

# Checks for some programs.
AC_PROG_INSTALL

# Checks for Fortran compilers.
FCFLAGS_PRESERVE=$FCFLAGS
AC_LANG(Fortran)
AC_PROG_FC(gfortran g95 ifort pgf95 xlf90_r xlf95_r)
FCFLAGS=$FCFLAGS_PRESERVE
# Restore FCFLAGS, because AC_PROG_FC() sets FCFLAGS to '-g -O2', '-O2' or '-g' undesirably.
# Desired FCFLAGS are:
case "$FC" in
  gfortran|g95)
    test -n "$FCFLAGS" || FCFLAGS="-g -Wall -O3 -funroll-loops -pipe -fopenmp"
    ;;
  ifort)
    test -n "$FCFLAGS" || FCFLAGS="-r8 -O3 -warn all -openmp -openmp-report2 -parallel -par-report1"
    ;;
  pgf*)
    test -n "$FCFLAGS" || FCFLAGS="-O3"
    ;;
  xlf*)
    test -n "$FCFLAGS" || FCFLAGS="-q64 -O4 -qsmp=omp -qinfo=all -qreport -qnoescape"
    LDFLAGS="$LDFLAGS -L/usr/local/lib -Wl,-blpdata -Wl,-bmaxdata:90000000000"
    ;;
  *)
    test -n "$FCFLAGS" || FCFLAGS="-O2"
    ;;
esac

## FFT library
AC_ARG_WITH(fft,[AS_HELP_STRING([--with-fft=ARG],[FFT library.
                                   ARG must be one of
                                   fftw3_omp (default, fast),
                                   fftw3_threads,
                                   mkl (Intel MKL, FFTW3 wrappers to Intel MKL functions),
                                   ssl2 (Fujitsu SSL II), or
                                   no (error).])])
case "$with_fft" in
  fftw3_omp | fftw | fftw3 | yes | "")
    AC_CHECK_LIB(fftw3, dfftw_execute,,
                 [AC_MSG_ERROR([You need FFTW3 library, i.e. libfftw3. See config.log.])])
    AC_CHECK_LIB(fftw3_omp, dfftw_init_threads,,
                 [AC_MSG_ERROR([You need FFTW3 library compiled with --enable-openmp, i.e. libfftw3_omp. See config.log.])])
    ;;
  fftw3_threads)
    AC_CHECK_LIB(fftw3, dfftw_execute,,
                 [AC_MSG_ERROR([You need FFTW3 library, i.e. libfftw3. See config.log.])])
    AC_CHECK_LIB(fftw3_threads, dfftw_init_threads,,
                 [AC_MSG_ERROR([You need FFTW3 library compiled with --enable-threads, i.e. libfftw3_threads. See config.log.])])
    ;;
  mkl)
    AC_CHECK_LIB(mkl_gf_lp64, dfftw_execute,[AC_DEFINE([HAVE_LIBFFTW3]) LIBS="-lmkl_core -lmkl_intel_thread -lmkl_gf_lp64 $LIBS"],
                 AC_MSG_ERROR([You need Intel Math Kernel Library (MKL). See config.log.]),[-lmkl_core -lmkl_intel_thread])
    ;;
  ssl2)
    ;;
  cufft*)
    AC_CHECK_LIB(libcufft, cufftPlan3d,, AC_MSG_ERROR([You need CUFFT. See config.log.]))
    ;;
  no)
    AC_MSG_ERROR([You need an FFT library.])
    ;;
  *)
    AC_MSG_ERROR([Unknown FFT library: $with_fft])
    ;;
esac

## LAPACK
AC_ARG_WITH(lapack,[AS_HELP_STRING([--with-lapack=ARG],[LAPACK library.
                                     ARG must be one of
                                     lapack (default),
                                     mkl (Intel MKL),
                                     ssl2 (Fujitsu SSL II).])])
case "$with_lapack" in
  lapack | yes | "")
    AC_CHECK_LIB(blas,  dsymv,,[AC_CHECK_LIB(f77blas,  dsymv,,AC_MSG_ERROR([You need  BLAS  library. See config.log.]))])
    AC_CHECK_LIB(lapack,dsyev,,AC_MSG_ERROR([You need LAPACK library. See config.log.]))
    ;;
  essl)
    AC_CHECK_LIB(esslsmp,dsyevx,,AC_MSG_ERROR([You need ESSL library for LAPACK subroutines. See config.log.]))
    ;;
  mkl)
    case "$LIBS" in
      *mkl_gf_lp64*)
        AC_MSG_NOTICE([lapack... It seems that there is -lmkl_gf_lp64 already in LIBS.])
        ;;
      *)
        AC_MSG_NOTICE([lapack... Add -lmkl_core -lmkl_intel_thread -lmkl_gf_lp64 at the tail of LIBS.])
        LIBS="$LIBS -lmkl_core -lmkl_intel_thread -lmkl_gf_lp64"
        ;;
    esac
    ;;
  ssl2)
    ;;
  *)
    AC_MSG_ERROR([Unknown LAPACK library: $with_lapack])
    ;;
esac

## free-form
case "$FC" in
  gfortran*)
    FCFLAGS="$FCFLAGS -ffree-form -ffree-line-length-none"
    ;;
  frt*)
    FCFLAGS="$FCFLAGS -Free"
    ;;
  xlf*)
    FCFLAGS="$FCFLAGS -qfree=f90"
    ;;
  *)
    AC_FC_FREEFORM()
    ;;
esac

# Finalize.
AC_SUBST(CPPFLAGS)

AC_CONFIG_FILES([Makefile
             src/Makefile])
AC_OUTPUT
