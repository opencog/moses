ADD_CXXTEST(vertexUTest)
TARGET_LINK_LIBRARIES(vertexUTest
	comboreduct
	comboant
	mosesutil
)

ADD_CXXTEST(messageUTest)
TARGET_LINK_LIBRARIES(messageUTest
	comboreduct
	comboant
	mosesutil
)

ADD_CXXTEST(listUTest)
TARGET_LINK_LIBRARIES(listUTest
	comboreduct
	comboant
	mosesutil
)

ADD_CXXTEST(similarityUTest)
TARGET_LINK_LIBRARIES(similarityUTest
	comboreduct
	comboant
	mosesutil
)

ADD_CXXTEST(type_treeUTest)
TARGET_LINK_LIBRARIES(type_treeUTest
	comboreduct
	comboant
	mosesutil
)

ADD_CXXTEST(annUTest)
TARGET_LINK_LIBRARIES(annUTest
	comboreduct
	mosesutil
)

ADD_CXXTEST(variable_unifierUTest)
TARGET_LINK_LIBRARIES(variable_unifierUTest
	comboreduct
	comboant
	mosesutil
)

ADD_CXXTEST(procedure_callUTest)
TARGET_LINK_LIBRARIES(procedure_callUTest
	comboreduct
	comboant
	mosesutil
)

ADD_CXXTEST(procedure_repositoryUTest)
TARGET_LINK_LIBRARIES(procedure_repositoryUTest
	comboreduct
	comboant
	mosesutil
)

ADD_CXXTEST(interpreterUTest)
TARGET_LINK_LIBRARIES(interpreterUTest
	comboreduct
	comboant
	mosesutil
)

ADD_CXXTEST(iostream_comboUTest)
TARGET_LINK_LIBRARIES(iostream_comboUTest
	comboreduct
	comboant
	mosesutil
)
