
ENABLE_TESTING()
INCLUDE(AddCxxtest)

IF (WIN32)
	ADD_DEFINITIONS(-D_CRT_SECURE_CPP_OVERLOAD_STANDARD_NAMES
	                -D_CRT_SECURE_NO_WARNINGS)
ENDIF (WIN32)

ADD_DEFINITIONS(-DPROJECT_SOURCE_DIR="${CMAKE_SOURCE_DIR}"
                -DPROJECT_BINARY_DIR="${CMAKE_BINARY_DIR}")

# Is this actually needed ???
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR})

# Perform tests in component-dependency order, as much as possible.
# For example, since everything depends on util and config functions,
# these are tested first. Since most things depend on the atomspace,
# that is tested next, etc.
IF (CXXTEST_FOUND)
	# ADD_SUBDIRECTORY (util)

	#ADD_SUBDIRECTORY (hillclimbing)

	ADD_SUBDIRECTORY (comboreduct)
	ADD_SUBDIRECTORY (moses)
	ADD_SUBDIRECTORY (feature-selection)

ENDIF (CXXTEST_FOUND)

IF (HAVE_CYTHON)
	ADD_SUBDIRECTORY (cython)
ENDIF (HAVE_CYTHON)
