#
# Master MOSES CMake file.
#
# General organization:
# -- check for different compilers, OS'es
# -- search for various required & optional libraries/tools
# -- decide what to build based on above results.
# -- configure various config files.
# -- print pretty summary
#

# MOSES' version, must be updated manually
SET (MOSES_VERSION_MAJOR 3)
SET (MOSES_VERSION_MINOR 6)
SET (MOSES_VERSION_PATCH 10)
SET (MOSES_VERSION "${MOSES_VERSION_MAJOR}.${MOSES_VERSION_MINOR}.${MOSES_VERSION_PATCH}")

# Need to have cmake-2.8.8 in order to be able to check the
# compiler version.
CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)
IF (COMMAND CMAKE_POLICY)
	CMAKE_POLICY(SET CMP0003 NEW)
	CMAKE_POLICY(SET CMP0005 OLD)
ENDIF (COMMAND CMAKE_POLICY)

IF(CMAKE_VERSION VERSION_GREATER 3.0.2)
	CMAKE_POLICY(SET CMP0037 OLD)   # In order to allow "make test"
ENDIF(CMAKE_VERSION VERSION_GREATER 3.0.2)

PROJECT(moses)

# ----------------------------------------------------------
# User-modifiable options. Feel free to change these!
#
# Uncomment to be in Release mode [default].
# SET(CMAKE_BUILD_TYPE Release)

# Uncomment to build in debug mode.
# SET(CMAKE_BUILD_TYPE Debug)

# Uncomment to be in coverage testing mode.
# SET(CMAKE_BUILD_TYPE Coverage)

# Uncomment to build in profile mode.
# SET(CMAKE_BUILD_TYPE Profile)

# uncomment to build in release mode with debug information
# SET(CMAKE_BUILD_TYPE RelWithDebInfo)

# Default build type
IF (CMAKE_BUILD_TYPE STREQUAL "")
	SET(CMAKE_BUILD_TYPE Release)
ENDIF (CMAKE_BUILD_TYPE STREQUAL "")

MESSAGE(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

ADD_DEFINITIONS(-DPROJECT_SOURCE_DIR=\\"${CMAKE_SOURCE_DIR}\\"
                -DPROJECT_BINARY_DIR=\\"${CMAKE_BINARY_DIR}\\")

# Add the 'lib' dir to cmake's module search path
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/lib/")

# ===============================================================
# Detect different compilers and OS'es, tweak flags as necessary.

IF (CMAKE_COMPILER_IS_GNUCXX)
	IF (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.6.4)
		MESSAGE(FATAL_ERROR "GCC version must be at least 4.6.4!")
	ENDIF (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.6.4)

	IF (APPLE)
		SET(CMAKE_C_FLAGS "-Wall -Wno-long-long -Wno-conversion")
		SET(CMAKE_C_FLAGS_DEBUG "-O0 -g")
		SET(CMAKE_C_FLAGS_PROFILE "-O0 -pg")
		SET(CMAKE_C_FLAGS_RELEASE "-O2 -g0")
		# Vital to do this otherwise unresolved symbols everywhere:
		SET(CMAKE_SHARED_LINKER_FLAGS "-Wl,-flat_namespace,-undefined,dynamic_lookup")
		SET(CMAKE_EXE_LINKER_FLAGS "-Wl,-flat_namespace,-undefined,dynamic_lookup")
	ELSE (APPLE)
		SET(CMAKE_C_FLAGS "-Wall -fPIC")
		SET(CMAKE_C_FLAGS "-Wl,--copy-dt-needed-entries")
		SET(CMAKE_C_FLAGS_DEBUG "-ggdb3 -fstack-protector")
		SET(CMAKE_C_FLAGS_PROFILE "-O2 -g3 -fstack-protector -pg")
		SET(CMAKE_C_FLAGS_RELEASE "-O2 -g -fstack-protector")
	ENDIF (APPLE)

	# 1) -Wno-variadic-macros is to avoid warnings regarding using
	# variadic in macro OC_ASSERT (the warning warns that this is only
	# available from C99, lol!)
	#
	# 2) -fopenmp for multithreading support
	#
	# 3) -std=gnu++0x for C++0x and GNU extensions support
	SET(CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS} -Wno-variadic-macros -fopenmp -std=gnu++0x")

	SET(CMAKE_CXX_FLAGS_DEBUG ${CMAKE_C_FLAGS_DEBUG})
	SET(CMAKE_CXX_FLAGS_PROFILE ${CMAKE_C_FLAGS_PROFILE})
	SET(CMAKE_CXX_FLAGS_RELEASE ${CMAKE_C_FLAGS_RELEASE})

	# Options for generating gcov code coverage output
	SET(CMAKE_C_FLAGS_COVERAGE "-O0 -g -fprofile-arcs -ftest-coverage -fno-inline")
	SET(CMAKE_CXX_FLAGS_COVERAGE "${CMAKE_C_FLAGS_COVERAGE} -fno-default-inline")
	# Might be needed for some combinations of ln and gcc
	IF (CMAKE_BUILD_TYPE STREQUAL "Coverage")
		LINK_LIBRARIES(gcov)
	ENDIF (CMAKE_BUILD_TYPE STREQUAL "Coverage")
ENDIF (CMAKE_COMPILER_IS_GNUCXX)

# Do the windows build
IF (WIN32)
	ADD_DEFINITIONS(-DWIN32)
	IF (CYGWIN)
		ADD_DEFINITIONS(-DCYGWIN)
	ENDIF (CYGWIN)

	# It is "not unix" when the code is compiled under windows but not
	# under cygwin.
	IF (NOT UNIX)
		ADD_DEFINITIONS(-DWIN32_NOT_UNIX)
		FIND_PACKAGE(PThreads REQUIRED)
		FIND_PACKAGE(STLPort REQUIRED)
		INCLUDE_DIRECTORIES("${CMAKE_SOURCE_DIR}/include/win32" ${PTHREADS_INCLUDE_DIR} ${STLPORT_INCLUDE_DIR})
		LINK_LIBRARIES(${PTHREADS_LIBRARY} ${STLPORT_LIBRARIES})
		ADD_DEFINITIONS(-D_CRT_SECURE_NO_WARNINGS -D_CRT_NONSTDC_NO_DEPRECATE)
	ELSE (NOT UNIX)
		LINK_LIBRARIES(pthread)
		IF (CYGWIN)
			SET(CMAKE_SHARED_LINKER_FLAGS "-Wl")
		ELSE (CYGWIN)
			SET(CMAKE_SHARED_LINKER_FLAGS "-Wl,--enable-new-dtags")
		ENDIF (CYGWIN)
	ENDIF (NOT UNIX)
ENDIF (WIN32)

# ===================================================================
# Check for existance of various required, optional packages.
# Listed in alphabetical order, more or less.

# Check for boost. We need dynamic-linked, threaded libs by default.
SET(Boost_USE_STATIC_LIBS OFF)
SET(Boost_USE_MULTITHREADED ON)

# Required boost packages
# hypertable requires boost iostreams
# boost-1.49 no longer has a libboost_iostreams
# 1.46 is minimum for required filesystem support
# program_options needed by some combo utilities
FIND_PACKAGE(Boost 1.46 COMPONENTS date_time filesystem program_options regex serialization system thread REQUIRED)

IF(Boost_FOUND)
	SET(Boost_FOUND_SAVE 1)

	# Boost.Signals2 requires this in signal_type.hpp
	SET(BOOST_PARAMETER_MAX_ARITY 7)
	ADD_DEFINITIONS(-DBOOST_PARAMETER_MAX_ARITY=${BOOST_PARAMETER_MAX_ARITY})
ELSE(Boost_FOUND)
	MESSAGE(FATAL_ERROR "Boost 1.46 or newer is needed to build OpenCog!")
ENDIF(Boost_FOUND)

MESSAGE(STATUS "Boost version ${Boost_VERSION} found.")

# Optional boost packages; can build without these.
FIND_PACKAGE(Boost 1.46 COMPONENTS python program_options QUIET)
FIND_PACKAGE(Boost 1.46 COMPONENTS math_c99 QUIET)

# Arghhh. Except cmake is treating above as required, not optional. #$%**&
IF(Boost_FOUND_SAVE)
	SET(Boost_FOUND 1)
ENDIF(Boost_FOUND_SAVE)

IF(Boost_PROGRAM_OPTIONS_FOUND)
	MESSAGE(STATUS "Found Boost::program_options")
ELSE(Boost_PROGRAM_OPTIONS_FOUND)
	MESSAGE(STATUS "Boost program_options missing: needed for Spatial tools (MapTools).")
ENDIF(Boost_PROGRAM_OPTIONS_FOUND)

IF(Boost_MATH_C99_FOUND)
	MESSAGE(STATUS "Found Boost::math")
ELSE(Boost_MATH_C99_FOUND)
	MESSAGE(STATUS "Boost math missing: needed to run combo learning tests.")
ENDIF(Boost_MATH_C99_FOUND)

# ----------------------------------------------------------
# Cogutil
FIND_PACKAGE(CogUtil CONFIG)
IF (COGUTIL_FOUND)
	MESSAGE(STATUS "CogUtil found.")
	ADD_DEFINITIONS(-DHAVE_COGUTIL)
	SET(HAVE_COGUTIL 1)
ELSE (COGUTIL_FOUND)
	MESSAGE(FATAL_ERROR "CogUtil missing: it is needed!")
ENDIF (COGUTIL_FOUND)

# add the 'cmake' directory from cogutil to search path
list(APPEND CMAKE_MODULE_PATH  ${COGUTIL_DATA_DIR}/cmake)

include(${COGUTIL_DATA_DIR}/cmake/Summary.cmake)

# ----------------------------------------------------------
FIND_PACKAGE(Cxxtest)
IF (CXXTEST_FOUND)
	MESSAGE(STATUS "Found CxxTest.")
ELSE (CXXTEST_FOUND)
	MESSAGE(STATUS "CxxTest missing: needed for unit tests.")
ENDIF (CXXTEST_FOUND)

# ----------------------------------------------------------
# MPI is needed for MOSES on compute clusters.
FIND_PACKAGE(MPI)

IF (MPI_FOUND)
	ADD_DEFINITIONS(-DHAVE_MPI)
	SET(HAVE_MPI 1)
	MESSAGE(STATUS "MPI was found.")
	INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
ELSE (MPI_FOUND)
	MESSAGE(STATUS "MPI missing: needed for MPI-MOSES.")
ENDIF (MPI_FOUND)

# ----------------------------------------------------------
# Python bindings using Cython
# Cython is used to generate bindings
FIND_PACKAGE(PythonLibs 2.7)
IF(PYTHONLIBS_FOUND)
	MESSAGE(STATUS "Python libs found.")
ELSE(PYTHONLIBS_FOUND)
	MESSAGE(STATUS "Python libs NOT found.")
ENDIF(PYTHONLIBS_FOUND)

FIND_PACKAGE(Cython 0.19.0)
IF (CYTHON_FOUND AND PYTHONLIBS_FOUND)
	# find python destination dir for python bindings
	ADD_DEFINITIONS(-DHAVE_CYTHON)
	SET(HAVE_CYTHON 1)
	EXECUTE_PROCESS(
		COMMAND ${PYTHON_EXECUTABLE} -c "from distutils.sysconfig import *;print(get_python_lib())"
		OUTPUT_VARIABLE PYTHON_DEST
	)

	# replace new line at end
	STRING(REPLACE "\n" "" PYTHON_DEST "${PYTHON_DEST}")
	IF ("${PYTHON_DEST}" STREQUAL "")
		MESSAGE(FATAL_ERROR "Python destination dir not found")
	ELSE ("${PYTHON_DEST}" STREQUAL "")
		MESSAGE(STATUS "Python destination dir found: ${PYTHON_DEST}" )
	ENDIF ("${PYTHON_DEST}" STREQUAL "")
ELSE (CYTHON_FOUND AND PYTHONLIBS_FOUND)
	IF(NOT CYTHON_FOUND)
		MESSAGE(STATUS "Cython executable not found.")
	ENDIF(NOT CYTHON_FOUND)
ENDIF (CYTHON_FOUND AND PYTHONLIBS_FOUND)

# Nosetests will find and automatically run python tests.
FIND_PROGRAM(NOSETESTS_EXECUTABLE nosetests)
FIND_PROGRAM(NOSETESTS3_EXECUTABLE nosetests3)
IF (CYTHON_FOUND AND PYTHONLIBS_FOUND)
	IF (NOSETESTS_EXECUTABLE)
		SET(HAVE_NOSETESTS 1)
	ENDIF (NOSETESTS_EXECUTABLE)

	IF (NOSETESTS3_EXECUTABLE)
		SET(HAVE_NOSETESTS 1)
		SET(NOSETESTS_EXECUTABLE ${NOSETESTS3_EXECUTABLE})
	ENDIF (NOSETESTS3_EXECUTABLE)
ENDIF (CYTHON_FOUND AND PYTHONLIBS_FOUND)

# ----------------------------------------------------------
# Optional, currently needed only to hush up DRD in util/Logger.cc
FIND_PACKAGE(VALGRIND)
IF (VALGRIND_FOUND)
	MESSAGE(STATUS "VALGRIND was found.")
	IF (VALGRIND_INCLUDE_DIR)
		MESSAGE(STATUS "VALGRIND devel headers found.")
		ADD_DEFINITIONS(-DHAVE_VALGRIND)
	ELSE (VALGRIND_INCLUDE_DIR)
		MESSAGE(STATUS "VALGRIND devel headers NOT FOUND: needed for thread debugging.")
	ENDIF (VALGRIND_INCLUDE_DIR)
ELSE (VALGRIND_FOUND)
	MESSAGE(STATUS "VALGRIND missing: needed for thread debugging.")
ENDIF (VALGRIND_FOUND)

# ==========================================================
# Decide what to build, based on the packages found.

IF(Boost_FOUND AND COGUTIL_FOUND)
	SET(HAVE_COMBOREDUCT 1)
	SET(HAVE_MOSES 1)
	SET(HAVE_FEATURE_SELECTION 1)
ENDIF(Boost_FOUND AND COGUTIL_FOUND)

# ===================================================================
# global includes

# set confdir and datadir
IF (NOT DEFINED CONFDIR)
	SET (CONFDIR "${CMAKE_INSTALL_PREFIX}/etc")
ENDIF (NOT DEFINED CONFDIR)
IF (NOT DEFINED DATADIR)
	SET (DATADIR "${CMAKE_INSTALL_PREFIX}/share/moses")
ENDIF (NOT DEFINED DATADIR)
ADD_DEFINITIONS(-DCONFDIR=\\"${CONFDIR}\\")
ADD_DEFINITIONS(-DDATADIR=\\"${DATADIR}\\")

# (re?)define MAN_INSTALL_DIR
SET (MAN_INSTALL_DIR "${CMAKE_INSTALL_PREFIX}/share/man")

# small hack to handle unixes that use "/usr/lib64" instead of "/usr/lib" as the
# default lib path on 64 bit archs
IF (NOT DEFINED LIB_DIR_SUFFIX)
	EXECUTE_PROCESS(COMMAND ${CMAKE_CXX_COMPILER} -print-search-dirs OUTPUT_VARIABLE PRINT_SEARCH_DIRS_OUTPUT)
	STRING(REGEX MATCH "\r?\nlibraries:.*\r?\n" COMPILER_LIB_SEARCH_DIRS ${PRINT_SEARCH_DIRS_OUTPUT})
	IF (NOT ${COMPILER_LIB_SEARCH_DIRS} STREQUAL "")
		STRING(REGEX MATCH "/lib64/:|/lib64:|/lib64\n" HAS_LIB64 ${COMPILER_LIB_SEARCH_DIRS})
		IF (NOT ${HAS_LIB64} STREQUAL "")
			SET(LIB_DIR_SUFFIX "64")
		ENDIF (NOT ${HAS_LIB64} STREQUAL "")
	ENDIF (NOT ${COMPILER_LIB_SEARCH_DIRS} STREQUAL "")
ENDIF (NOT DEFINED LIB_DIR_SUFFIX)

# set default include paths
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR} ${Boost_INCLUDE_DIRS}
	${COGUTIL_INCLUDE_DIR})

# rpath handling
SET(CMAKE_SKIP_BUILD_RPATH	FALSE)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib/moses")
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

ADD_SUBDIRECTORY(moses)

IF (CXXTEST_FOUND AND COGUTIL_FOUND)
	ADD_CUSTOM_TARGET(tests)
	ADD_SUBDIRECTORY(tests EXCLUDE_FROM_ALL)
	IF (CMAKE_BUILD_TYPE STREQUAL "Coverage")
		# doing coverage stuff while running tests if this is the Coverage build
		ADD_CUSTOM_TARGET(test
			# TODO lcov should be found by cmake first
			# TODO set it up so that we can pick to run coverage per test, or
			# combined across all tests (the latter is MUCH faster). Use a define?
			# There is coverage specific stuff in AddCxxTest.cmake now...
			# -
			# Depends on cogserver because RESTFulTest needs it
			# and cmake has no way to add non-test dependencies to tests
			DEPENDS tests cogserver
			WORKING_DIRECTORY tests
			COMMAND ${CMAKE_CTEST_COMMAND} --force-new-ctest-process $(ARGS)
			# This script combines the coverage analysis of each test,
			# then creates html in tests/lcov
			# Note: this should now be run separately...
			#COMMAND ${PROJECT_SOURCE_DIR}/scripts/combine_lcov.sh
			COMMENT "Running tests with coverage..."
		)
	ELSE (CMAKE_BUILD_TYPE STREQUAL "Coverage")
		# If this is a build with coverage enabled then test normally
		ADD_CUSTOM_TARGET(test
			DEPENDS tests
			WORKING_DIRECTORY tests
			COMMAND ${CMAKE_CTEST_COMMAND} --force-new-ctest-process $(ARGS)
			COMMENT "Running tests..."
		)
	ENDIF (CMAKE_BUILD_TYPE STREQUAL "Coverage")
ENDIF (CXXTEST_FOUND AND COGUTIL_FOUND)


ADD_SUBDIRECTORY(examples EXCLUDE_FROM_ALL)

IF (NOT WIN32)
	ADD_CUSTOM_TARGET (examples
		# using CMAKE_BUILD_TOOL results in teh cryptic error message
		# warning: jobserver unavailable: using -j1.  Add `+' to parent make rule.
		# This is because make doesn't know how to pass jobserver args to
		# the submake.  So, instead, just use $(MAKE) (with round parens)
		# -- that will do the right thing.
		# COMMAND ${CMAKE_BUILD_TOOL}
		COMMAND $(MAKE)
		WORKING_DIRECTORY examples
		COMMENT "Building examples"
	)
ENDIF (NOT WIN32)


ADD_CUSTOM_TARGET(cscope
	COMMAND find moses tests -name '*.cc' -o -name '*.h' -o -name '*.cxxtest' > ${CMAKE_SOURCE_DIR}/cscope.files
	COMMAND cscope -b
	WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
	COMMENT "Generating CScope database"
)

# packaging
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Meta-Optimizing Semantic Evolutionary Search")
SET(CPACK_PACKAGE_NAME "moses")
SET(CPACK_PACKAGE_VENDOR "opencog.org")
SET(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/README.md")
SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE")
SET(CPACK_PACKAGE_VERSION_MAJOR "0")
SET(CPACK_PACKAGE_VERSION_MINOR "1")
SET(CPACK_PACKAGE_VERSION_PATCH "4")
SET(CPACK_BINARY_STGZ "OFF")
SET(CPACK_BINARY_TBZ2 "OFF")
SET(CPACK_BINARY_TGZ "OFF")
SET(CPACK_BINARY_TZ "OFF")
SET(CPACK_SOURCE_STGZ "OFF")
SET(CPACK_SOURCE_TBZ2 "ON")
SET(CPACK_SOURCE_TGZ "OFF")
SET(CPACK_SOURCE_TZ "OFF")
SET(CPACK_GENERATOR "TBZ2")
SET(CPACK_SOURCE_GENERATOR "TBZ2")
SET(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")
SET(CPACK_SOURCE_PACKAGE_FILE_NAME "${CPACK_PACKAGE_FILE_NAME}")
SET(CPACK_SOURCE_IGNORE_FILES
	"~$"
	".a$"
	".so$"
	".log$"
	".tar.z$"
	".tar.gz$"
	".tar.bz2$"
	"/Makefile$"
	"/CMakeFiles/"
	"/CMakeCache.txt$"
	"/DartTestfile.txt$"
	"/cmake_install.cmake$"
	"/Testfile.cmake$"
	"/CTestTestfile.cmake$"
	"/CTestTestfile.cmake$"
	"UTest$"
	"UTest.cpp$"
	"/_CPack_Packages/"
	"/CPackConfig.cmake$"
	"/CPackSourceConfig.cmake$"
	"/sniff$"
	"/install_manifest.txt$"
	"${CPACK_PACKAGE_FILE_NAME}"
	"^${PROJECT_SOURCE_DIR}/debian/"
	"^${PROJECT_SOURCE_DIR}/.*spec$"
	"^${PROJECT_SOURCE_DIR}/vendor/"
	"^${PROJECT_SOURCE_DIR}/.bzr/"
	"^${PROJECT_SOURCE_DIR}/.bzrignore$"
	"^${PROJECT_SOURCE_DIR}/.git/"
)
SET(CPACK_PACKAGE_INSTALL_DIRECTORY "CMake ${CMake_VERSION_MAJOR}.${CMake_VERSION_MINOR}")
IF (WIN32)
	# Its "not unix" is cygwin is being used.
	IF (NOT UNIX)
		# There is a bug in NSI that does not handle full unix paths
		# properly. Make sure there is at least one set of four (4)
		# backslashes.
		SET(CPACK_PACKAGE_ICON "${CMake_SOURCE_DIR}\\\\lib\\\\opencog.ico")
		SET(CPACK_NSIS_DISPLAY_NAME "${CPACK_PACKAGE_INSTALL_DIRECTORY} The Open Cognition Framework")
		SET(CPACK_NSIS_HELP_LINK "http:\\\\\\\\www.opencog.org")
		SET(CPACK_NSIS_URL_INFO_ABOUT "http:\\\\\\\\www.opencog.org")
		SET(CPACK_NSIS_CONTACT "opencog@googlegroups.com")
		SET(CPACK_NSIS_MODIFY_PATH ON)
	ELSE(NOT UNIX)
		SET(CPACK_STRIP_FILES "bin/moses")
		SET(CPACK_SOURCE_STRIP_FILES "")
	ENDIF(NOT UNIX)
ENDIF(WIN32)

SET(CPACK_PACKAGE_EXECUTABLES "cogserver" "The Open Cognition Framework")
INCLUDE(CPack)

# documentation
FIND_PACKAGE(Doxygen)
ADD_SUBDIRECTORY(doc EXCLUDE_FROM_ALL)
ADD_SUBDIRECTORY(lib)

# Show a summary of what we got
SUMMARY_ADD("ComboReduct" "Library for reduction of combo program trees" HAVE_COMBOREDUCT)
SUMMARY_ADD("Doxygen" "Code documentation" DOXYGEN_FOUND)
SUMMARY_ADD("Feature selection" "MOSES based feature selection tool" HAVE_FEATURE_SELECTION)
SUMMARY_ADD("MOSES" "Meta-Optimizing Semantic Evolutionary Search" HAVE_MOSES)
SUMMARY_ADD("Python bindings" "Python (cython) bindings" HAVE_CYTHON)
SUMMARY_ADD("Python tests" "Python bindings nose tests" HAVE_NOSETESTS)
SUMMARY_ADD("Unit tests" "Unit tests" CXXTEST_FOUND AND COGUTIL_FOUND)
SUMMARY_SHOW()
