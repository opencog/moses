
ADD_EXECUTABLE (pleasure-main main.cpp)
TARGET_LINK_LIBRARIES (pleasure-main
	pleasure
	comboreduct
	comboant
	${COGUTIL_LIBRARY}
)

ADD_CUSTOM_COMMAND(
	TARGET main
	POST_BUILD
COMMAND ${CMAKE_COMMAND} -E copy ${OPENCOG_SOURCE_DIR}/moses/learning/pleasure/main/nodelist ${OPENCOG_BINARY_DIR}/moses/learning/pleasure/main
)


