#
# Tests are listed and run in more-or-less dependency order:
# Later tests assume earlier tests pass.
#
ADD_CXXTEST(iostream_scored_combo_treeUTest)
TARGET_LINK_LIBRARIES(iostream_scored_combo_treeUTest
	moses
	comboreduct
	${COGUTIL_LIBRARY}
)

ADD_CXXTEST(KnobTurningUTest)
TARGET_LINK_LIBRARIES(KnobTurningUTest
	moses
	comboreduct
	${COGUTIL_LIBRARY}
)

# requires boost::math; don't build if not found.
IF(Boost_MATH_C99_FOUND)
	ADD_CXXTEST(NeighborSamplingUTest)
	TARGET_LINK_LIBRARIES(NeighborSamplingUTest
		moses
		comboreduct
		${COGUTIL_LIBRARY}
	)
ENDIF(Boost_MATH_C99_FOUND)


ADD_CXXTEST(FieldSetUTest)
TARGET_LINK_LIBRARIES(FieldSetUTest
	moses
	comboreduct
	${COGUTIL_LIBRARY}
)

#ADD_CXXTEST(representationUTest)
#TARGET_LINK_LIBRARIES(representationUTest
#	moses
#	comboreduct
#	${COGUTIL_LIBRARY}
#)

ADD_CXXTEST(scoringUTest)
TARGET_LINK_LIBRARIES(scoringUTest
	moses
	comboreduct
	${COGUTIL_LIBRARY}
)

ADD_CXXTEST(metapopulationUTest)
TARGET_LINK_LIBRARIES(metapopulationUTest
	moses
	comboreduct
	${COGUTIL_LIBRARY}
)

ADD_CXXTEST(optimizationUTest)
TARGET_LINK_LIBRARIES(optimizationUTest
	moses
	comboreduct
	${COGUTIL_LIBRARY}
)

ADD_CXXTEST(MOSESUTest)
TARGET_LINK_LIBRARIES(MOSESUTest
	moses_exec
	moses
	comboreduct
	${COGUTIL_LIBRARY}
)

#ADD_CXXTEST(weightedUTest)
#TARGET_LINK_LIBRARIES(weightedUTest
#	moses_exec
#	moses
#	comboreduct
#	${COGUTIL_LIBRARY}
#)

#ADD_CXXTEST(mixedUTest)
#TARGET_LINK_LIBRARIES(mixedUTest
#	moses_exec
#	moses
#	comboreduct
#	${COGUTIL_LIBRARY}
#)

#ADD_CXXTEST(irisUTest)
#TARGET_LINK_LIBRARIES(irisUTest
#	moses_exec
#	moses
#	comboreduct
#	${COGUTIL_LIBRARY}
#)

#ADD_CXXTEST(selectionUTest)
#TARGET_LINK_LIBRARIES(selectionUTest
#	moses_exec
#	moses
#	comboreduct
#	${COGUTIL_LIBRARY}
#)

#ADD_CXXTEST(dynFeatSelUTest)
#TARGET_LINK_LIBRARIES(dynFeatSelUTest
#	moses_exec
#	moses
#	comboreduct
#	${COGUTIL_LIBRARY}
#)

ADD_CXXTEST(AntUTest)
TARGET_LINK_LIBRARIES(AntUTest
	ant_scoring
	moses_exec
	moses
	comboreduct
	${COGUTIL_LIBRARY}
)

#ADD_CXXTEST(diversityUTest)
#TARGET_LINK_LIBRARIES(diversityUTest
#	moses_exec
#	moses
#	comboreduct
#	${COGUTIL_LIBRARY}
#)
