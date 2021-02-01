cmake_minimum_required(VERSION 3.9)

project(Mesh)

set (CMAKE_EXPORT_COMPILE_COMMANDS ON)
if( EXISTS "${CMAKE_CURRENT_BINARY_DIR}/compile_commands.json" )
  EXECUTE_PROCESS( COMMAND ${CMAKE_COMMAND} -E copy_if_different
    ${CMAKE_CURRENT_BINARY_DIR}/compile_commands.json
    ${CMAKE_CURRENT_SOURCE_DIR}/compile_commands.json
  )
endif()

if(NOT CMAKE_BUILD_TYPE)
    set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
            "Choose the type of build, options are: None Debug Release." FORCE)
endif(NOT CMAKE_BUILD_TYPE)

option(BUILD_PYTHON_INTERFACE "Build the python interface" ON)

add_subdirectory(cpp)
if(BUILD_PYTHON_INTERFACE)
   add_subdirectory(py)
endif()
