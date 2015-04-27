#
# Tests are listed and run in more-or-less dependency order:
# Later tests assume earlier tests pass.
#
ADD_CXXTEST(iostream_scored_combo_treeUTest)
TARGET_LINK_LIBRARIES(iostream_scored_combo_treeUTest
	moses
	comboreduct
	mosesutil
)

ADD_CXXTEST(KnobTurningUTest)
TARGET_LINK_LIBRARIES(KnobTurningUTest
	moses
	comboreduct
	mosesutil
)

# requires boost::math; don't build if not found.
IF(Boost_MATH_C99_FOUND)
	ADD_CXXTEST(NeighborSamplingUTest)
	TARGET_LINK_LIBRARIES(NeighborSamplingUTest
		moses
		comboreduct
		mosesutil
	)
ENDIF(Boost_MATH_C99_FOUND)


ADD_CXXTEST(FieldSetUTest)
TARGET_LINK_LIBRARIES(FieldSetUTest
	moses
	comboreduct
	mosesutil
)

ADD_CXXTEST(representationUTest)
TARGET_LINK_LIBRARIES(representationUTest
	moses
	comboreduct
	mosesutil
)

ADD_CXXTEST(scoringUTest)
TARGET_LINK_LIBRARIES(scoringUTest
	moses
	comboreduct
	mosesutil
)

ADD_CXXTEST(metapopulationUTest)
TARGET_LINK_LIBRARIES(metapopulationUTest
	moses
	comboreduct
	mosesutil
)

ADD_CXXTEST(optimizationUTest)
TARGET_LINK_LIBRARIES(optimizationUTest
	moses
	comboreduct
	mosesutil
)

ADD_CXXTEST(MOSESUTest)
TARGET_LINK_LIBRARIES(MOSESUTest
	moses_exec
	moses
	comboreduct
	mosesutil
)

ADD_CXXTEST(weightedUTest)
TARGET_LINK_LIBRARIES(weightedUTest
	moses_exec
	moses
	comboreduct
	mosesutil
)

ADD_CXXTEST(mixedUTest)
TARGET_LINK_LIBRARIES(mixedUTest
	moses_exec
	moses
	comboreduct
	mosesutil
)

ADD_CXXTEST(irisUTest)
TARGET_LINK_LIBRARIES(irisUTest
	moses_exec
	moses
	comboreduct
	mosesutil
)

ADD_CXXTEST(selectionUTest)
TARGET_LINK_LIBRARIES(selectionUTest
	moses_exec
	moses
	comboreduct
	mosesutil
)

ADD_CXXTEST(dynFeatSelUTest)
TARGET_LINK_LIBRARIES(dynFeatSelUTest
	moses_exec
	moses
	comboreduct
	mosesutil
)

ADD_CXXTEST(AntUTest)
TARGET_LINK_LIBRARIES(AntUTest
	ant_scoring
	moses_exec
	moses
	comboreduct
	mosesutil
)
