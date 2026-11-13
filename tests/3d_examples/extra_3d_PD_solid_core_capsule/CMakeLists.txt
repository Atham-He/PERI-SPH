cmake_minimum_required(VERSION 3.16)
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    project(PDSolidCoreCapsule LANGUAGES CXX)
    set(SPHINXSYS_2D OFF CACHE BOOL "" FORCE)
    set(SPHINXSYS_3D ON CACHE BOOL "" FORCE)
    set(SPHINXSYS_BUILD_TESTS OFF CACHE BOOL "" FORCE)
    set(SPHINXSYS_USE_ONEDPL OFF CACHE BOOL "" FORCE)
    set(SPHINXSYS_DEVELOPER_MODE OFF CACHE BOOL "" FORCE)
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../../.. sphinxsys EXCLUDE_FROM_ALL)
endif()
set(PERI_SHARED ${CMAKE_CURRENT_SOURCE_DIR}/../../../src/shared)
file(SHA256 ${CMAKE_CURRENT_SOURCE_DIR}/pd_solid_capsule.cpp MAIN_HASH)
file(SHA256 ${CMAKE_CURRENT_SOURCE_DIR}/capsule_core_fixture.h FIXTURE_HASH)
file(SHA256 ${PERI_SHARED}/peridynamics/pd_solid_reference_coupling.h COUPLING_HASH)
file(SHA256 ${PERI_SHARED}/materials/drucker_prager_cap_kernel.h KERNEL_HASH)
file(SHA256 ${PERI_SHARED}/materials/drucker_prager_cap_solid.h ADAPTER_HASH)
set(PD_ASSEMBLY_HASHES "${MAIN_HASH};${FIXTURE_HASH};${COUPLING_HASH}")
foreach(PD_HEADER pd_quadratic_correspondence.h pd_principal_stress_onset.h
        pd_fracture_correspondence.h pd_fracture_surface.h pd_fragment_contact.h pd_energy_crack_band.h
        pd_layered_geometry.h pd_layered_correspondence.h pd_layered_surface.h pd_coulomb_friction.h)
    set(PD_HEADER_PATH ${PERI_SHARED}/peridynamics/${PD_HEADER})
    set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS ${PD_HEADER_PATH})
    file(SHA256 ${PD_HEADER_PATH} PD_HEADER_HASH)
    string(APPEND PD_ASSEMBLY_HASHES ";${PD_HEADER_HASH}")
endforeach()
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS
    ${CMAKE_CURRENT_SOURCE_DIR}/pd_solid_capsule.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/capsule_core_fixture.h
    ${PERI_SHARED}/peridynamics/pd_solid_reference_coupling.h
    ${PERI_SHARED}/materials/drucker_prager_cap_kernel.h
    ${PERI_SHARED}/materials/drucker_prager_cap_solid.h)
string(SHA256 ASSEMBLY_HASH "${PD_ASSEMBLY_HASHES}")
string(SHA256 MATERIAL_HASH "${KERNEL_HASH};${ADAPTER_HASH}")
add_executable(pd_solid_capsule pd_solid_capsule.cpp)
target_link_libraries(pd_solid_capsule sphinxsys_3d)
target_compile_features(pd_solid_capsule PRIVATE cxx_std_17)
target_compile_definitions(pd_solid_capsule PRIVATE
    CAPSULE_SOURCE_SHA256="${ASSEMBLY_HASH}"
    DPC_MATERIAL_SHA256="${MATERIAL_HASH}")

add_executable(pd_coulomb_law_check test_coulomb_friction.cpp)
target_link_libraries(pd_coulomb_law_check sphinxsys_3d)
target_compile_features(pd_coulomb_law_check PRIVATE cxx_std_17)

add_executable(dpc_moving_junction_check test_dpc_moving_junction.cpp)
target_link_libraries(dpc_moving_junction_check sphinxsys_3d)
target_compile_features(dpc_moving_junction_check PRIVATE cxx_std_17)
