cmake_minimum_required(VERSION 3.16)
# A dedicated build avoids configuring or running unrelated FSI studies.
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    project(PDCantileverVerification LANGUAGES CXX)
    set(SPHINXSYS_2D OFF CACHE BOOL "" FORCE)
    set(SPHINXSYS_3D ON CACHE BOOL "" FORCE)
    set(SPHINXSYS_BUILD_TESTS OFF CACHE BOOL "" FORCE)
    set(SPHINXSYS_USE_ONEDPL OFF CACHE BOOL "" FORCE)
    set(SPHINXSYS_DEVELOPER_MODE OFF CACHE BOOL "" FORCE)
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../../.. sphinxsys EXCLUDE_FROM_ALL)
    enable_testing()
endif()
STRING(REGEX REPLACE ".*/(.*)" "\\1" CURRENT_FOLDER ${CMAKE_CURRENT_SOURCE_DIR})
PROJECT("${CURRENT_FOLDER}")

SET(EXECUTABLE_OUTPUT_PATH "${PROJECT_BINARY_DIR}/bin/")
SET(BUILD_INPUT_PATH "${EXECUTABLE_OUTPUT_PATH}/input")
SET(BUILD_RELOAD_PATH "${EXECUTABLE_OUTPUT_PATH}/reload")

file(MAKE_DIRECTORY ${BUILD_INPUT_PATH})

add_executable(${PROJECT_NAME})
target_sources(${PROJECT_NAME} PRIVATE PD_cantilever.cpp)
target_link_libraries(${PROJECT_NAME} sphinxsys_3d)
set_target_properties(${PROJECT_NAME} PROPERTIES VS_DEBUGGER_WORKING_DIRECTORY "${EXECUTABLE_OUTPUT_PATH}")

add_executable(pd_cantilever_benchmark cantilever_benchmark.cpp)
target_link_libraries(pd_cantilever_benchmark sphinxsys_3d)
file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/smoke)
add_test(NAME pd_cantilever_smoke
    COMMAND pd_cantilever_benchmark --case-config=${CMAKE_CURRENT_SOURCE_DIR}/smoke.cfg --state_recording=false
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/smoke)
set_tests_properties(pd_cantilever_smoke PROPERTIES TIMEOUT 120 LABELS "pd;smoke")
