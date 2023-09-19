# We need to put this in a library, so that AntUTest can link to it.
ADD_LIBRARY (ant_scoring
	ant_scoring.cc
)

ADD_EXECUTABLE (moses-ant-hillclimbing
	moses-ant-hillclimbing.cc
	ant_scoring.cc
)

ADD_EXECUTABLE (moses-ant-particleswarm
	moses-ant-particleswarm.cc
	ant_scoring.cc
)

TARGET_LINK_LIBRARIES (moses-ant-hillclimbing
	moses
	comboreduct
	comboant
	${COGUTIL_LIBRARY}
)

TARGET_LINK_LIBRARIES (moses-ant-particleswarm
	moses
	comboreduct
	comboant
	${COGUTIL_LIBRARY}
)
