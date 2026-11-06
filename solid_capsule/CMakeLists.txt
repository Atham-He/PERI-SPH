cmake_minimum_required(VERSION 3.18)
project(SolidCapsuleSPH LANGUAGES CXX)
get_filename_component(SPHINXSYS_PROJECT_DIR "${CMAKE_CURRENT_LIST_DIR}/.." ABSOLUTE)
set(SPHINXSYS_2D OFF CACHE BOOL "" FORCE)
set(SPHINXSYS_3D ON CACHE BOOL "" FORCE)
set(SPHINXSYS_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(SPHINXSYS_USE_ONEDPL OFF CACHE BOOL "" FORCE)
set(SPHINXSYS_DEVELOPER_MODE OFF CACHE BOOL "" FORCE)
add_subdirectory("${SPHINXSYS_PROJECT_DIR}" sphinxsys EXCLUDE_FROM_ALL)
# GCC 13 on this workstation ICEs in fmt's fwprop1 in this existing CK source.
# Keep the workaround local to this build entry; do not change the library source.
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 13
   AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 14)
    set_property(SOURCE "${SPHINXSYS_PROJECT_DIR}/src/shared/shared_ck/particle_dynamics/solid_dynamics/derived_solid_state.cpp"
        DIRECTORY "${SPHINXSYS_PROJECT_DIR}/src" APPEND PROPERTY COMPILE_OPTIONS "-O0")
endif()
enable_testing()
add_subdirectory("${SPHINXSYS_PROJECT_DIR}/tests/3d_examples/test_3d_solid_core_shell_compression" capsule)
