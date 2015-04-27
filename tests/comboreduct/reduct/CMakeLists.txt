#ADD_CXXTEST(ActionReductUTest)
#TARGET_LINK_LIBRARIES(ActionReductUTest
#	comboreduct
#   util
#)

ADD_CXXTEST(BooleanReductUTest)
TARGET_LINK_LIBRARIES(BooleanReductUTest
	comboreduct
	comboant
	mosesutil
)

ADD_CXXTEST(BooleanFullUTest)
TARGET_LINK_LIBRARIES(BooleanFullUTest
	comboreduct
	comboant
	mosesutil
)

ADD_CXXTEST(ContinReductUTest)
TARGET_LINK_LIBRARIES(ContinReductUTest
	comboreduct
	comboant
	mosesutil
)

ADD_CXXTEST(FullReductUTest)
TARGET_LINK_LIBRARIES(FullReductUTest
	comboreduct
	comboant
	mosesutil
)

ADD_CXXTEST(MixedReductUTest)
TARGET_LINK_LIBRARIES(MixedReductUTest
	comboreduct
	comboant
	mosesutil
)

#ADD_CXXTEST(PerceptionReductUTest)
#TARGET_LINK_LIBRARIES(PerceptionReductUTest
#	comboreduct
#   util
#)

ADD_CXXTEST(ANNReductUTest)
TARGET_LINK_LIBRARIES(ANNReductUTest
	comboreduct
	comboant
	mosesutil
)

ADD_CXXTEST(FoldReductUTest)
TARGET_LINK_LIBRARIES(FoldReductUTest
	comboreduct
	comboant
	mosesutil
)
