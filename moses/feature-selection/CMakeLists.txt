ADD_SUBDIRECTORY(scorers)
ADD_SUBDIRECTORY(main)
ADD_SUBDIRECTORY(man)

ADD_LIBRARY(feature_selection SHARED
	algo/deme_optimize.cc
	algo/incremental.cc
	algo/random.cc
	algo/simple.cc
	algo/stochastic_max_dependency.cc
	scorers/moses_optim.cc
	main/feature-selection.cc
)

TARGET_LINK_LIBRARIES(feature_selection
	# INTERFACE moses
	${COGUTIL_LIBRARY}
	${Boost_PROGRAM_OPTIONS_LIBRARY} 
)

# Install library
IF (WIN32)
	INSTALL(TARGETS feature_selection
	        EXPORT MosesTargets
	        DESTINATION "lib${LIB_DIR_SUFFIX}/moses")
ELSE (WIN32)
	INSTALL(TARGETS feature_selection
	        EXPORT MosesTargets
	        LIBRARY DESTINATION "lib${LIB_DIR_SUFFIX}"    # lib*.so files
	        ARCHIVE DESTINATION "lib${LIB_DIR_SUFFIX}")   # lib*.a files
ENDIF (WIN32)
