add_library(Math INTERFACE)
add_library(Mesh::Math ALIAS Math)

target_compile_features(Math INTERFACE cxx_std_17)
target_compile_definitions(Math INTERFACE Math=1)
target_include_directories(Math
  INTERFACE
        $<INSTALL_INTERFACE:include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        )
#set_target_properties(Math PROPERTIES OUTPUT_NAME "math$<$<CONFIG:Debug>:-d>")

set(Math_Version 1.0)

include(GNUInstallDirs)
install(TARGETS Math
        EXPORT MathTargets
        LIBRARY  DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE  DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME  DESTINATION ${CMAKE_INSTALL_BINDIR}
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(EXPORT MathTargets
        FILE MathTargets.cmake
        NAMESPACE Mesh::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Mesh)

include(CMakePackageConfigHelpers)
write_basic_package_version_file("MathConfigVersion.cmake"
        VERSION ${Math_Version}
        COMPATIBILITY SameMajorVersion)

configure_package_config_file(${CMAKE_CURRENT_LIST_DIR}/cmake/MathConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/MathConfig.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Mesh)

install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/MathConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/MathConfigVersion.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Mesh)
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

export(TARGETS Math FILE ${CMAKE_CURRENT_BINARY_DIR}/MathTargets.cmake NAMESPACE Mesh::)
export(PACKAGE Math)
