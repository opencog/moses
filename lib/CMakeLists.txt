CONFIGURE_FILE(${CMAKE_SOURCE_DIR}/lib/moses.conf
    ${PROJECT_BINARY_DIR}/lib/moses.conf)
# Create a test configuration file by copying the default configuration file
CONFIGURE_FILE(${CMAKE_SOURCE_DIR}/lib/moses-test.conf
    ${PROJECT_BINARY_DIR}/lib/moses-test.conf)

export(EXPORT MosesTargets
          FILE "${CMAKE_CURRENT_BINARY_DIR}/Moses/MosesTargets.cmake"
      )

set(ConfigPackageLocation lib/cmake/Moses)

install(EXPORT MosesTargets
  FILE
    MosesTargets.cmake
  DESTINATION
    ${ConfigPackageLocation}
)

include(CMakePackageConfigHelpers)
configure_package_config_file(MosesConfig.cmake.in ${CMAKE_CURRENT_BINARY_DIR}/MOSESConfig.cmake
                              INSTALL_DESTINATION lib/Moses/cmake
                              PATH_VARS CMAKE_INSTALL_PREFIX)


write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/MOSESConfigVersion.cmake"
  VERSION ${MOSES_VERSION}
  COMPATIBILITY SameMajorVersion
)

INSTALL(FILES
       ${CMAKE_CURRENT_BINARY_DIR}/MOSESConfigVersion.cmake
       ${CMAKE_CURRENT_BINARY_DIR}/MOSESConfig.cmake
       DESTINATION
       ${ConfigPackageLocation})
