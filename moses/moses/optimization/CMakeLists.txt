#install header files
INSTALL(FILES
	hill-climbing.h
	optimization.h
	star-anneal.h
	univariate.h
	particle-swarm.h
	hybrid-hc-ps.h

	DESTINATION

	"include/${PROJECT_NAME}/moses/optimization" 
)
