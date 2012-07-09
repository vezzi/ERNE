dnl Process this file with autoconf to produce a configure script.

AC_PREREQ([2.67])
AC_INIT([ERNE],[1.0.1])

AC_LANG(C++)

AC_CONFIG_MACRO_DIR([m4])

AC_CANONICAL_TARGET
AM_INIT_AUTOMAKE()
LT_INIT

# standard system requirements
AC_PROG_CC
AC_PROG_CXX
AM_PROG_CC_C_O
AC_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_SIZE_T
AC_CHECK_HEADERS([pthread.h stdio.h math.h limits.h bzlib.h zlib.h iostream fstream string vector exception sstream], [], [AC_MSG_ERROR(["Header not found"])] )
AC_PROG_INSTALL
AC_PROG_MAKE_SET

# common flags
CXXFLAGS+=" -Wall"
CFLAGS+=" -Wall -D_USE_KNETFILE -D_FILE_OFFSET_BITS=64"

# checking for 64 bit system (optional)
echo $ECHO_N "checking for 64 bit architecture... "
case "$build_cpu" in
x86_64*)
	CXXFLAGS+=" -m64 -DCPU_64"
	CFLAGS+=" -m64 -DCPU_64"
	COMPILING_64="yes"
	echo "yes"
	;;
*) echo "no" ;;
esac

# checking for BOOST library (required)
echo "===== BOOST support (required) ====="
BOOST_STATIC()
BOOST_REQUIRE([1.40.0])
BOOST_IOSTREAMS
BOOST_PROGRAM_OPTIONS
BOOST_THREADS
BOOST_REGEX

# enable debug programs
AC_ARG_ENABLE([maintainer-mode],
     [AS_HELP_STRING([--enable-maintainer-mode],[create debug programs [no]], [28])],
     [enable_debug=$enableval],
     [enable_debug=no])
AM_CONDITIONAL(DEBUG_ENABLED,test "$enable_debug" != "no")

# enable static programs
AC_ARG_ENABLE([static-programs],
     [AS_HELP_STRING([--enable-static-programs],[create static program version [no]], [28])],
     [enable_static_programs=$enableval],
     [enable_static_programs=no])
AM_CONDITIONAL(STATIC_ENABLED,test "$enable_static_programs" != "no")

# check for MPI (optional)
AC_ARG_ENABLE([mpi],
     [AS_HELP_STRING([--enable-mpi],[create OpenMPI based program(s) [check]])],
     [enable_mpi=$enableval],
     [enable_mpi=no]) # default: check

if test "$enable_mpi" != "no" ; then
	echo "===== OpenMPI support (optional) ====="
	ACX_MPI()
	if test "$MPICXX" != "g++" ; then
		AC_LANG_PUSH(C)
		ACX_MPI()
		AC_LANG_POP()
		if test "$MPICC" != "cpp" ; then
			MPI=yes
		fi
	fi
fi

AM_CONDITIONAL(HAVE_MPI,test "$MPI" = "yes")
AM_CONDITIONAL(CYGWIN_TARGET,test "$target_os" = "cygwin" )

# print summary
echo     "==================== Summary ====================="
if test "$COMPILING_64" = "yes" ; then
	echo "= Compile with 64 bit                            =";
else
	echo "= Compile with 32 bit                            =";
fi
if test "$MPI" = "yes"; then
	echo "= Compile with MPI                               =";
else
	echo "= Compile without MPI                            =";
fi
if test "$enable_debug" != "no" ; then
	echo "= Compile with debug programs                    =";
else
	echo "= Compile without debug programs                 =";
fi
if test "$enable_static_programs" != "no" ; then
	echo "= Compile static programs version (no ERNE-DMAP) =";
else
	echo "= Do not compile static programs version         =";
fi
if test "$target_os" = "cygwin" ; then
	echo "= Detected WINDOWS system                        =";
fi
echo     "=================================================="

# write files
AC_CONFIG_FILES(Makefile src/Makefile)
AC_OUTPUT

