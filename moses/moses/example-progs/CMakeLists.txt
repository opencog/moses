ADD_EXECUTABLE (moses-ann-pole1 moses-ann-pole1)
TARGET_LINK_LIBRARIES (moses-ann-pole1
	moses
	mosesutil
)

ADD_EXECUTABLE (moses-ann-pole2 moses-ann-pole2 pole_balancing)
TARGET_LINK_LIBRARIES (moses-ann-pole2
	moses
	mosesutil
)

ADD_EXECUTABLE (moses-ann-pole2nv moses-ann-pole2nv pole_balancing)
TARGET_LINK_LIBRARIES (moses-ann-pole2nv
	moses
	mosesutil
)

ADD_EXECUTABLE (moses-ann-pole2-hillclimbing
	moses-ann-pole2-hillclimbing
	pole_balancing)
TARGET_LINK_LIBRARIES (moses-ann-pole2-hillclimbing
	moses
	mosesutil
)

ADD_EXECUTABLE (moses-ann-xor moses-ann-xor)
TARGET_LINK_LIBRARIES (moses-ann-xor
	moses
	mosesutil
)

ADD_EXECUTABLE (continmax continmax)
TARGET_LINK_LIBRARIES (continmax
	moses
	mosesutil
)

ADD_EXECUTABLE (nmax nmax)
TARGET_LINK_LIBRARIES (nmax
	moses
	mosesutil
)

ADD_EXECUTABLE (onemax onemax)
TARGET_LINK_LIBRARIES (onemax
	moses
	mosesutil
)

ADD_EXECUTABLE (ontomax ontomax)
TARGET_LINK_LIBRARIES (ontomax
	moses
	mosesutil
)

ADD_EXECUTABLE (trap-bit trap-bit)
TARGET_LINK_LIBRARIES (trap-bit
	moses
	mosesutil
)

ADD_EXECUTABLE (trap-uni trap-uni)
TARGET_LINK_LIBRARIES (trap-uni
	moses
	mosesutil
)

ADD_EXECUTABLE (moses-perf moses-perf)
TARGET_LINK_LIBRARIES (moses-perf
   moses_exec
	moses
	mosesutil
)
