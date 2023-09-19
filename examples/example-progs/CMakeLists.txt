ADD_EXECUTABLE (moses-ann-pole1 moses-ann-pole1.cc)
TARGET_LINK_LIBRARIES (moses-ann-pole1
	moses
	${COGUTIL_LIBRARY}
)

ADD_EXECUTABLE (moses-ann-pole2 moses-ann-pole2.cc pole_balancing.cc)
TARGET_LINK_LIBRARIES (moses-ann-pole2
	moses
	${COGUTIL_LIBRARY}
)

ADD_EXECUTABLE (moses-ann-pole2nv moses-ann-pole2nv.cc pole_balancing.cc)
TARGET_LINK_LIBRARIES (moses-ann-pole2nv
	moses
	${COGUTIL_LIBRARY}
)

ADD_EXECUTABLE (moses-ann-pole2-hillclimbing
	moses-ann-pole2-hillclimbing.cc
	pole_balancing.cc)
TARGET_LINK_LIBRARIES (moses-ann-pole2-hillclimbing
	moses
	${COGUTIL_LIBRARY}
)

ADD_EXECUTABLE (moses-ann-xor moses-ann-xor.cc)
TARGET_LINK_LIBRARIES (moses-ann-xor
	moses
	${COGUTIL_LIBRARY}
)

ADD_EXECUTABLE (continmax continmax.cc)
TARGET_LINK_LIBRARIES (continmax
	moses
	${COGUTIL_LIBRARY}
)

ADD_EXECUTABLE (nmax nmax.cc)
TARGET_LINK_LIBRARIES (nmax
	moses
	${COGUTIL_LIBRARY}
)

ADD_EXECUTABLE (onemax onemax.cc)
TARGET_LINK_LIBRARIES (onemax
	moses
	${COGUTIL_LIBRARY}
)

ADD_EXECUTABLE (ontomax ontomax.cc)
TARGET_LINK_LIBRARIES (ontomax
	moses
	${COGUTIL_LIBRARY}
)

ADD_EXECUTABLE (trap-bit trap-bit.cc)
TARGET_LINK_LIBRARIES (trap-bit
	moses
	${COGUTIL_LIBRARY}
)

ADD_EXECUTABLE (trap-uni trap-uni.cc)
TARGET_LINK_LIBRARIES (trap-uni
	moses
	${COGUTIL_LIBRARY}
)

ADD_EXECUTABLE (moses-perf moses-perf.cc)
TARGET_LINK_LIBRARIES (moses-perf
   moses_exec
	moses
	${COGUTIL_LIBRARY}
)
