add_executable(extra_3d_PD_gel_compression gel_compression.cpp)
target_link_libraries(extra_3d_PD_gel_compression PRIVATE sphinxsys_3d SPHinXsys::peri_execution)
set_target_properties(extra_3d_PD_gel_compression PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/bin")

# This test parses/validates inputs only: no particles, output or solver launch.
add_test(NAME peri_gel_model_description COMMAND extra_3d_PD_gel_compression
    --case-config "${CMAKE_CURRENT_SOURCE_DIR}/configs/numerical_control.cfg" --describe-model)
set_tests_properties(peri_gel_model_description PROPERTIES LABELS peri_fast TIMEOUT 10
    PASS_REGULAR_EXPRESSION "solid_sphere_no_separate_liquid")

find_package(Python3 COMPONENTS Interpreter REQUIRED)
add_test(NAME peri_gel_preflight COMMAND ${Python3_EXECUTABLE}
    "${SPHINXSYS_PROJECT_DIR}/tools/peri_gel_preflight.py"
    --executable "$<TARGET_FILE:extra_3d_PD_gel_compression>")
set_tests_properties(peri_gel_preflight PROPERTIES LABELS peri_fast TIMEOUT 30)
