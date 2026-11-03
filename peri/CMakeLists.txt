cmake_minimum_required(VERSION 3.16)
project(PeriIntegration LANGUAGES CXX)

# A separate build entry configures only the retained Peri examples.
get_filename_component(SPHINXSYS_PROJECT_DIR "${CMAKE_CURRENT_LIST_DIR}/.." ABSOLUTE)
set(SPHINXSYS_2D OFF CACHE BOOL "" FORCE)
set(SPHINXSYS_3D ON CACHE BOOL "" FORCE)
set(SPHINXSYS_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(SPHINXSYS_USE_ONEDPL OFF CACHE BOOL "" FORCE)
set(SPHINXSYS_DEVELOPER_MODE OFF CACHE BOOL "" FORCE)
add_subdirectory("${SPHINXSYS_PROJECT_DIR}" sphinxsys EXCLUDE_FROM_ALL)
enable_testing()
find_package(GTest CONFIG REQUIRED)
option(PERI_ENABLE_LONG_TESTS "Enable retained long Peri application tests" OFF)

foreach(example IN ITEMS
    extra_3d_PD_cantilever
    extra_3d_PD_capsule_compression
    extra_3d_PD_gel_compression
    extra_3d_PD_framework_tests
    extra_3d_PD_brittle_disc
    extra_3d_PD_ice_jet)
    add_subdirectory("${SPHINXSYS_PROJECT_DIR}/tests/3d_examples/${example}" "${example}")
endforeach()
