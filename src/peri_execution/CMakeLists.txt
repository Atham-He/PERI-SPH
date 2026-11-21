set(PERI_EXECUTION_HEADERS "${CMAKE_CURRENT_SOURCE_DIR}/../shared/peridynamics")

add_library(sphinxsys_peri_execution STATIC pd_execution.cpp)
add_library(SPHinXsys::peri_execution ALIAS sphinxsys_peri_execution)
target_compile_features(sphinxsys_peri_execution PUBLIC cxx_std_17)
target_include_directories(sphinxsys_peri_execution PUBLIC
    $<BUILD_INTERFACE:${PERI_EXECUTION_HEADERS}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../shared/materials>)

if(PERI_SPH_ENABLE_CUDA)
    find_package(CUDAToolkit REQUIRED)
    target_sources(sphinxsys_peri_execution PRIVATE
        "${PERI_EXECUTION_HEADERS}/pd_crack_projection_cuda.cu"
        "${PERI_EXECUTION_HEADERS}/pd_bonded_interface_cuda.cu"
        "${PERI_EXECUTION_HEADERS}/pd_spatial_derivatives_cuda.cu"
        dpc_point_cuda.cu
        pd_solid_time_step_cuda.cu
        pd_nosb_deformation_cuda.cu
        pd_fluid_traction_cuda.cu
        pd_fluid_acoustic_cuda.cu
        pd_fluid_viscosity_cuda.cu
        pd_fluid_cfl_cuda.cu)
    set_target_properties(sphinxsys_peri_execution PROPERTIES CUDA_STANDARD 17 CUDA_STANDARD_REQUIRED YES)
    target_compile_options(sphinxsys_peri_execution PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:--fmad=false>)
    target_compile_definitions(sphinxsys_peri_execution PRIVATE PERI_SPH_HAS_CUDA=1)
    set_source_files_properties(dpc_point_cuda.cu "${PERI_EXECUTION_HEADERS}/pd_spatial_derivatives_cuda.cu"
        PROPERTIES COMPILE_OPTIONS "--expt-relaxed-constexpr")
    target_link_libraries(sphinxsys_peri_execution PRIVATE CUDA::cudart_static Eigen3::Eigen)
endif()
