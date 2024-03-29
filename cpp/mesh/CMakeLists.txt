add_library(Mesh SHARED mesh.cc segment.cc elements.cc system.cc meshing.cc)
add_library(Mesh::Mesh ALIAS Mesh)

find_package(Eigen3 REQUIRED)
find_package(Math REQUIRED)
find_package(Triangle REQUIRED)
find_package(IO REQUIRED)

target_compile_features(Mesh PRIVATE cxx_std_17)
target_include_directories(Mesh
  PUBLIC
        $<INSTALL_INTERFACE:include>
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        PRIVATE
    ${EIGEN3_INCLUDE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}
)

get_target_property(Triangle_COMPILE_DEFINITIONS Mesh::Triangle COMPILE_DEFINITIONS)
target_compile_definitions(Mesh PRIVATE Triangle_COMPILE_DEFINITIONS)

set_target_properties(Mesh PROPERTIES OUTPUT_NAME "mesh$<$<CONFIG:Debug>:-d>")
#set_target_properties(Mesh PROPERTIES LINKER_LANGUAGE C)

set(Mesh_Version 1.0)
target_link_libraries(Mesh PRIVATE Mesh::IO Mesh::Math Mesh::Triangle)

include(GNUInstallDirs)
install(TARGETS Mesh
	EXPORT MeshTargets
        LIBRARY  DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE  DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME  DESTINATION ${CMAKE_INSTALL_BINDIR}
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(EXPORT MeshTargets
	FILE MeshTargets.cmake
        NAMESPACE Mesh::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Mesh)

include(CMakePackageConfigHelpers)
write_basic_package_version_file("MeshConfigVersion.cmake"
  VERSION ${Mesh_Version}
        COMPATIBILITY SameMajorVersion)

configure_package_config_file(${CMAKE_CURRENT_LIST_DIR}/cmake/MeshConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/MeshConfig.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Mesh)

install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/MeshConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/MeshConfigVersion.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Mesh)
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

export(TARGETS Mesh FILE ${CMAKE_CURRENT_BINARY_DIR}/MeshTargets.cmake NAMESPACE Mesh::)
export(PACKAGE Mesh)
