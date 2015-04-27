#install header files
INSTALL(FILES

    fs_scorer_base.h
    mutual_info.h
    moses_optim.h
    moses_matrix.h

	DESTINATION

	"include/${PROJECT_NAME}/feature-selection/scorers" 
)
