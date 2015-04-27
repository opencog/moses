
# We need to put this in a library, so that MOSESUTest can link to it.
ADD_LIBRARY (moses_exec
    demo-problems
    problem
    problem-params
    moses_exec
    moses_exec_def
    table-problems
)

ADD_EXECUTABLE (moses-exec moses-exec)
TARGET_LINK_LIBRARIES (moses-exec
	moses_exec
	moses
	comboant
)
SET_TARGET_PROPERTIES(moses-exec PROPERTIES OUTPUT_NAME "moses")

ADD_EXECUTABLE (eval-diversity eval-diversity)
TARGET_LINK_LIBRARIES (eval-diversity
	moses
	comboant
	mosesutil
	${Boost_PROGRAM_OPTIONS_LIBRARY} 
)

ADD_EXECUTABLE (eval-candidate eval-candidate)
TARGET_LINK_LIBRARIES (eval-candidate
	moses
	comboant
	mosesutil
	${Boost_PROGRAM_OPTIONS_LIBRARY}
)

ADD_EXECUTABLE (eval-candidate-likelihood eval-candidate-likelihood)
TARGET_LINK_LIBRARIES (eval-candidate-likelihood
	moses
	comboant
	mosesutil
	${Boost_PROGRAM_OPTIONS_LIBRARY}
)

# install moses-exec
INSTALL(TARGETS moses-exec eval-diversity eval-candidate eval-candidate-likelihood RUNTIME DESTINATION bin)

# Install library
IF (WIN32)
	INSTALL(TARGETS moses DESTINATION "lib${LIB_DIR_SUFFIX}/moses")
ELSE (WIN32)
	INSTALL(TARGETS moses_exec
	        LIBRARY DESTINATION "lib${LIB_DIR_SUFFIX}"    # lib*.so files
	        ARCHIVE DESTINATION "lib${LIB_DIR_SUFFIX}")   # lib*.a files
ENDIF (WIN32)

#install header files
INSTALL(FILES
	moses_exec.h
	moses_exec_def.h

	DESTINATION

	"include/${PROJECT_NAME}/moses/main" 
)
