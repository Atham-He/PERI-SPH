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
set(CRACK_CONTACT_IMPL ${PERI_SHARED}/../for_3D_build/particle_dynamics/solid_dynamics/pd_isoparametric_crack_surface.cpp)
file(SHA256 ${CRACK_CONTACT_IMPL} CRACK_CONTACT_HASH)
file(SHA256 ${PERI_SHARED}/materials/drucker_prager_cap_kernel.h KERNEL_HASH)
file(SHA256 ${PERI_SHARED}/materials/dpc_point_batch.h DPC_BATCH_HASH)
file(SHA256 ${PERI_SHARED}/materials/drucker_prager_cap_solid.h ADAPTER_HASH)
set(PD_ASSEMBLY_HASHES "${MAIN_HASH};${FIXTURE_HASH};${COUPLING_HASH};${CRACK_CONTACT_HASH}")
foreach(PD_HEADER pd_quadratic_correspondence.h pd_principal_stress_onset.h
        pd_fracture_correspondence.h pd_fracture_surface.h pd_fragment_contact.h pd_energy_crack_band.h pd_directional_crack_energy.h
        pd_layered_geometry.h pd_layered_correspondence.h pd_layered_surface.h pd_coulomb_friction.h
        pd_layered_fracture_reference.h pd_broken_bond_energy.h pd_fracture_visible_surface.h
        pd_cohesive_crack_law.h pd_cohesive_interface_imex.h pd_cohesive_surface_evolution.h
        pd_extrinsic_crack_growth.h pd_fracture_patch_partition.h pd_consistent_reconstruction.h pd_isoparametric_crack_surface.h pd_crack_face_contact.h pd_fracture_topology_checkpoint.h)
    set(PD_HEADER_PATH ${PERI_SHARED}/peridynamics/${PD_HEADER})
    set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS ${PD_HEADER_PATH})
    file(SHA256 ${PD_HEADER_PATH} PD_HEADER_HASH)
    string(APPEND PD_ASSEMBLY_HASHES ";${PD_HEADER_HASH}")
endforeach()
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS
    ${CMAKE_CURRENT_SOURCE_DIR}/pd_solid_capsule.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/capsule_core_fixture.h
    ${PERI_SHARED}/peridynamics/pd_solid_reference_coupling.h
    ${CRACK_CONTACT_IMPL}
    ${PERI_SHARED}/materials/drucker_prager_cap_kernel.h
    ${PERI_SHARED}/materials/dpc_point_batch.h
    ${PERI_SHARED}/materials/drucker_prager_cap_solid.h)
string(SHA256 ASSEMBLY_HASH "${PD_ASSEMBLY_HASHES}")
string(SHA256 MATERIAL_HASH "${KERNEL_HASH};${ADAPTER_HASH};${DPC_BATCH_HASH}")
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
target_link_libraries(dpc_moving_junction_check sphinxsys_3d SPHinXsys::peri_execution)
target_compile_features(dpc_moving_junction_check PRIVATE cxx_std_17)

# Optional S-NSNI-based bare-solid entry; default PD/SPH entry is unchanged.
add_executable(snsni_bare snsni_bare.cpp)
target_link_libraries(snsni_bare sphinxsys_3d SPHinXsys::peri_execution)
target_compile_features(snsni_bare PRIVATE cxx_std_17)

# S-NSNI solid core coupled to the layered NOSB-PD capsule shell.
set(SNSNI_CAPSULE_HASHES "${ASSEMBLY_HASH};${MATERIAL_HASH}")
foreach(SNSNI_FILE
        ${CMAKE_CURRENT_SOURCE_DIR}/pd_snsni_capsule.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/pd_snsni_capsule_fixture.h
        ${CMAKE_CURRENT_SOURCE_DIR}/recovery_checkpoint.h
        ${PERI_SHARED}/particle_dynamics/solid_dynamics/solid_snsni_plastic.h
        ${PERI_SHARED}/materials/dpc_snsni_adapter.h
        ${PERI_SHARED}/peridynamics/pd_radial_layered_surface.h
        ${PERI_SHARED}/peridynamics/pd_column_face_surface.h
        ${PERI_SHARED}/peridynamics/pd_execution.h
        ${PERI_SHARED}/../peri_execution/pd_execution.cpp
        ${PERI_SHARED}/peridynamics/pd_solid_time_step.h
        ${PERI_SHARED}/../peri_execution/pd_solid_time_step_cuda.cu
        ${PERI_SHARED}/../peri_execution/dpc_point_cuda.cu
        ${PERI_SHARED}/peridynamics/pd_spatial_derivatives.h
        ${PERI_SHARED}/peridynamics/pd_smoothed_force_projection.h
        ${PERI_SHARED}/peridynamics/pd_spatial_derivatives_cuda.cu
        ${PERI_SHARED}/peridynamics/pd_particle_trace_map.h
        ${PERI_SHARED}/peridynamics/pd_crack_projection.h
        ${PERI_SHARED}/peridynamics/pd_crack_projection_cuda.cu
        ${PERI_SHARED}/peridynamics/pd_bonded_interface_batch.h
        ${PERI_SHARED}/peridynamics/pd_bonded_interface_cuda.cu
        ${PERI_SHARED}/peridynamics/pd_cuda_buffer.cuh
        ${PERI_SHARED}/peridynamics/diffuse_crack_release.h)
    set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS ${SNSNI_FILE})
    file(SHA256 ${SNSNI_FILE} SNSNI_FILE_HASH)
    string(APPEND SNSNI_CAPSULE_HASHES ";${SNSNI_FILE_HASH}")
endforeach()
string(SHA256 SNSNI_CAPSULE_HASH "${SNSNI_CAPSULE_HASHES}")
add_executable(pd_snsni_capsule pd_snsni_capsule.cpp)
target_link_libraries(pd_snsni_capsule sphinxsys_3d SPHinXsys::peri_execution)
target_compile_features(pd_snsni_capsule PRIVATE cxx_std_17)
target_compile_definitions(pd_snsni_capsule PRIVATE
    CAPSULE_SOURCE_SHA256="${SNSNI_CAPSULE_HASH}"
    DPC_MATERIAL_SHA256="${MATERIAL_HASH}")

# Resolved NOSB-PD/solid-SPH capsule inside an explicit rigid-grain shear box.
file(SHA256 ${CMAKE_CURRENT_SOURCE_DIR}/pd_snsni_direct_shear.cpp DIRECT_SHEAR_MAIN_HASH)
file(SHA256 ${PERI_SHARED}/peridynamics/rigid_grain_shear.h DIRECT_SHEAR_CONTACT_HASH)
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS
    ${CMAKE_CURRENT_SOURCE_DIR}/pd_snsni_direct_shear.cpp
    ${PERI_SHARED}/peridynamics/rigid_grain_shear.h)
string(SHA256 DIRECT_SHEAR_HASH "${ASSEMBLY_HASH};${DIRECT_SHEAR_MAIN_HASH};${DIRECT_SHEAR_CONTACT_HASH}")
add_executable(pd_snsni_direct_shear pd_snsni_direct_shear.cpp)
target_link_libraries(pd_snsni_direct_shear sphinxsys_3d SPHinXsys::peri_execution)
target_compile_features(pd_snsni_direct_shear PRIVATE cxx_std_17)
target_compile_definitions(pd_snsni_direct_shear PRIVATE
    CAPSULE_SOURCE_SHA256="${DIRECT_SHEAR_HASH}"
    DPC_MATERIAL_SHA256="${MATERIAL_HASH}")

# Single frozen crack-layer acceptance: opening, closure, slip, exact transpose
# forces, physical explicit integration, and unchanged Gc work budget.
add_executable(pd_directional_crack_benchmark pd_directional_crack_benchmark.cpp)
target_link_libraries(pd_directional_crack_benchmark sphinxsys_3d)
target_compile_features(pd_directional_crack_benchmark PRIVATE cxx_std_17)
include(CTest)
if(BUILD_TESTING)
    add_test(NAME pd_directional_crack_open_close_slip
        COMMAND pd_directional_crack_benchmark ${CMAKE_CURRENT_BINARY_DIR}/directional_crack_acceptance)
    set_tests_properties(pd_directional_crack_open_close_slip PROPERTIES RUN_SERIAL TRUE)
endif()

# Actual disconnected precrack with one-sided reconstruction and face contact.
add_executable(pd_physical_crack_benchmark pd_physical_crack_benchmark.cpp)
target_link_libraries(pd_physical_crack_benchmark sphinxsys_3d SPHinXsys::peri_execution)
target_compile_features(pd_physical_crack_benchmark PRIVATE cxx_std_17)
if(BUILD_TESTING)
    add_test(NAME pd_physical_precrack_open_close_slip COMMAND pd_physical_crack_benchmark ${CMAKE_CURRENT_BINARY_DIR}/physical_precrack_acceptance)
    set_tests_properties(pd_physical_precrack_open_close_slip PROPERTIES RUN_SERIAL TRUE)
endif()
