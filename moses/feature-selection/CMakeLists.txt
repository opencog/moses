FIND_LIBRARY( MOSES_LIB moses  )

ADD_SUBDIRECTORY(scorers)
ADD_SUBDIRECTORY(main)
ADD_SUBDIRECTORY(man)

ADD_LIBRARY(feature_selection STATIC
	algo/deme_optimize
	algo/incremental
	algo/random
	algo/simple
	algo/stochastic_max_dependency
	scorers/moses_optim
	main/feature-selection
)

TARGET_LINK_LIBRARIES(feature_selection
    moses
	${COGUTIL_LIBRARY}
	${Boost_PROGRAM_OPTIONS_LIBRARY} 
)

# Install library
IF (WIN32)
	INSTALL(TARGETS feature_selection DESTINATION "lib${LIB_DIR_SUFFIX}/moses")
ELSE (WIN32)
	INSTALL(TARGETS feature_selection
	        LIBRARY DESTINATION "lib${LIB_DIR_SUFFIX}"    # lib*.so files
	        ARCHIVE DESTINATION "lib${LIB_DIR_SUFFIX}")   # lib*.a files
ENDIF (WIN32)
