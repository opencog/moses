# Need to use -fno-strict-aliasing when compiling cython code, in order
# to avoid nasty compiler warnings about aliasing.  Cython explicitly
# performs aliasing, in order to emulate python object inheritance.
# See, for example,
# https://groups.google.com/forum/#!topic/cython-users/JV1-KvIUeIg
#
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-strict-aliasing")

INCLUDE_DIRECTORIES(
	${PYTHON_INCLUDE_DIRS}
	${CMAKE_CURRENT_SOURCE_DIR}
	${CMAKE_CURRENT_BINARY_DIR}
)

SET(CYTHON_FLAGS "-f")

# Cythonizes one .pyx file into a .cpp file
# Additional arguments are dependencies
MACRO(CYTHON_ADD_MODULE_PYX name)
	IF(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${name}.pxd)
		SET(DEPENDS ${name}.pyx ${name}.pxd)
	ENDIF(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${name}.pxd)
	# Allow the user to specify dependencies as optional arguments
	SET(DEPENDS ${DEPENDS} ${ARGN})
	ADD_CUSTOM_COMMAND(
		OUTPUT ${name}.cpp
		COMMAND ${CYTHON_EXECUTABLE}
		ARGS ${CYTHON_FLAGS} -I ${PROJECT_BINARY_DIR}
		      -I ${ATOMSPACE_INCLUDE_DIR}/opencog/cython
		      # -I /usr/local/include/opencog/cython
		      -o ${name}.cpp
		      --cplus ${CMAKE_CURRENT_SOURCE_DIR}/${name}.pyx
		DEPENDS ${DEPENDS}
		COMMENT "Cythonizing ${name}.pyx")
ENDMACRO(CYTHON_ADD_MODULE_PYX)


# Use this as a guide:
# https://github.com/OpenKinect/libfreenect/blob/master/wrappers/python/CMakeLists.txt
IF (HAVE_MOSES)
	CYTHON_ADD_MODULE_PYX(pymoses
		# "<moses/moses/main/moses_exec.h>"
	)

	# opencog.pymoses Python bindings
	ADD_LIBRARY(moses_cython SHARED
		pymoses.cpp
	)

	# Don't reverse the order of moses_exec and moses, if reversed will give
	# undefined symbol: _ZNK5boost15program_options22error_with_option_name4whatEv
	# XXX Why???
	TARGET_LINK_LIBRARIES(moses_cython
		moses_exec
		moses
	)

	IF (HAVE_MPI)
		TARGET_LINK_LIBRARIES (moses_cython
			# These two work for both MPICH2 and OpenMPI, for me.
			${MPI_LIBRARY}
			${MPI_EXTRA_LIBRARY}
		)
	ENDIF (HAVE_MPI)

	SET_TARGET_PROPERTIES(moses_cython PROPERTIES
		PREFIX ""
		OUTPUT_NAME pymoses)

	INSTALL (TARGETS moses_cython
		EXPORT MosesTargets
		DESTINATION "${DATADIR}/python/moses")
ENDIF (HAVE_MOSES)

INSTALL (FILES
	__init__.py
DESTINATION "${DATADIR}/python/moses")
