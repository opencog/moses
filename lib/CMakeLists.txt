CONFIGURE_FILE(${CMAKE_SOURCE_DIR}/lib/moses.conf
    ${PROJECT_BINARY_DIR}/lib/moses.conf)
# Create a test configuration file by copying the default configuration file
CONFIGURE_FILE(${CMAKE_SOURCE_DIR}/lib/moses-test.conf
    ${PROJECT_BINARY_DIR}/lib/moses-test.conf)
