ADD_DEFINITIONS(-DMOSES_VERSION_MAJOR="${MOSES_VERSION_MAJOR}")
ADD_DEFINITIONS(-DMOSES_VERSION_MINOR="${MOSES_VERSION_MINOR}")
ADD_DEFINITIONS(-DMOSES_VERSION_PATCH="${MOSES_VERSION_PATCH}")

# MOSES git describe, automatically updated
EXECUTE_PROCESS(
	COMMAND git describe
	OUTPUT_VARIABLE MOSES_GIT_DESCRIBE
	WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

IF (MOSES_GIT_DESCRIBE)
  ADD_DEFINITIONS(-DMOSES_GIT_DESCRIBE="${MOSES_GIT_DESCRIBE}")
ENDIF (MOSES_GIT_DESCRIBE)

ADD_SUBDIRECTORY (moses)
ADD_SUBDIRECTORY (scoring)
ADD_SUBDIRECTORY (metapopulation)
ADD_SUBDIRECTORY (deme)
ADD_SUBDIRECTORY (representation)
ADD_SUBDIRECTORY (optimization)
ADD_SUBDIRECTORY (eda)
ADD_SUBDIRECTORY (main)
ADD_SUBDIRECTORY (man)

ADD_LIBRARY (moses SHARED
	deme/deme_expander.cc
	deme/feature_selector.cc

	eda/initialization.cc
	eda/local_structure.cc

	metapopulation/diversity.cc
	metapopulation/domination.cc
	metapopulation/ensemble.cc
	metapopulation/merging.cc
	metapopulation/metapopulation.cc
	metapopulation/metapop_params.cc
	metapopulation/subsample.cc

	moses/complexity.cc
	moses/distributed_moses.cc
	moses/local_moses.cc
	moses/moses_main.cc
	moses/mpi_moses.cc
	moses/neighborhood_sampling.cc
	moses/partial.cc
	moses/types.cc

	optimization/hill-climbing.cc
	optimization/optimization.cc
	optimization/star-anneal.cc
	optimization/univariate.cc
	optimization/particle-swarm.cc

	representation/build_knobs.cc
	representation/field_set.cc
	representation/knob_mapper.h
	representation/knobs.cc
	representation/representation.cc

	scoring/behave_cscore.cc
	scoring/bscores.cc
	scoring/discriminating_bscore.cc
	scoring/precision_bscore.cc
	scoring/bscores.cc
	scoring/time_dispersion.cc
	scoring/ss_bscore.cc
	scoring/scoring_base.cc
	scoring/select_bscore.cc
)

TARGET_LINK_LIBRARIES (moses
	feature_selection
	comboreduct
	comboant
	${COGUTIL_LIBRARY}
	${Boost_PROGRAM_OPTIONS_LIBRARY}
	${Boost_REGEX_LIBRARY}
)

IF (HAVE_MPI)
	TARGET_LINK_LIBRARIES (moses
		# These two work for both MPICH2 and OpenMPI, for me.
		${MPI_LIBRARY}
		${MPI_EXTRA_LIBRARY}
	)
ENDIF (HAVE_MPI)


# target_link_libraries(hello ${MPI_LIBRARIES})
# if(MPI_COMPILE_FLAGS)
#   set_target_properties(hello PROPERTIES
#     COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
# endif()
#
# if(MPI_LINK_FLAGS)
#   set_target_properties(hello PROPERTIES
#     LINK_FLAGS "${MPI_LINK_FLAGS}")
# endif()

# Install library
IF (WIN32)
	INSTALL(TARGETS moses
                EXPORT MosesTargets
                DESTINATION "lib${LIB_DIR_SUFFIX}/moses")
ELSE (WIN32)
	INSTALL(TARGETS moses
                EXPORT MosesTargets
	        LIBRARY DESTINATION "lib${LIB_DIR_SUFFIX}"    # lib*.so files
	        ARCHIVE DESTINATION "lib${LIB_DIR_SUFFIX}")   # lib*.a files
ENDIF (WIN32)
