#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.65])
AC_INIT([indri-pruning-toolkit], [1.0], [rueycheng@gmail.com])
AM_INIT_AUTOMAKE([foreign -Wall -Werror])
AC_CONFIG_SRCDIR([src/IndriPruneIndex.cpp])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CXX

# Checks for indri settings.
AC_ARG_WITH([indri],
	    [AS_HELP_STRING([--with-indri=PATH], [set the indri installation path [/usr/local]])],
	    [], [with_indri=/usr/local])
if test ! -d $with_indri; then
    AC_MSG_ERROR([Invalid path '$with_indri' (--with-indri)])
fi

indri_lib="$with_indri/lib"
indri_include="$with_indri/include"
if test ! -d $indri_lib; then 
    AC_MSG_ERROR([Path '$indri_lib' does not exist]); 
fi
if test ! -d $indri_include; then 
    AC_MSG_ERROR([Path '$indri_include' does not exist]); 
fi

LDFLAGS="$LDFLAGS -L$indri_lib"
CPPFLAGS="$CPPFLAGS -I$indri_include"
AC_SUBST([LDFLAGS])
AC_SUBST([CPPFLAGS])

# Checks for libraries.
AC_CHECK_LIB([m], [log])
AC_CHECK_LIB([pthread], [pthread_create])
AC_CHECK_LIB([z], [deflate])
AC_CHECK_LIB([iberty], [mkstemps])

AC_LANG_PUSH([C++])
AC_CHECK_LIB([indri], [main],
	     [], [AC_MSG_ERROR([cannot find the required library 'indri'])])
AC_LANG_POP([C++])

# Checks for header files.
AC_CHECK_HEADERS([indri/indri-platform.h],
		 [],[AC_MSG_ERROR([cannot find the required Indri headers])])

# Checks for typedefs, structures, and compiler characteristics.
AX_CXX_COMPILE_STDCXX_0X

if test "$ax_cv_cxx_compile_cxx0x_cxx" = yes; then
    CXXFLAGS="$CXXFLAGS -std=c++0x"
elif test "$ax_cv_cxx_compile_cxx0x_gxx" = yes; then 
    CXXFLAGS="$CXXFLAGS -std=gnu++0x"
elif test "$ax_cv_cxx_compile_cxx0x_native" = no; then
    AC_MSG_ERROR([cannot proceed without c++0x support])
fi

AX_BOOST_BASE([1.46],[],
	      [AC_MSG_ERROR([cannot find the boost libraries])])
AX_BOOST_FILESYSTEM
AX_BOOST_UNIT_TEST_FRAMEWORK

AC_SUBST([CXXFLAGS])

# Checks for library functions.

AC_CONFIG_FILES([Makefile src/Makefile test/Makefile])
AC_OUTPUT
