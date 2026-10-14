STRING(REGEX REPLACE ".*/(.*)" "\\1" CURRENT_FOLDER ${CMAKE_CURRENT_SOURCE_DIR})
PROJECT("${CURRENT_FOLDER}")

SET(EXECUTABLE_OUTPUT_PATH "${PROJECT_BINARY_DIR}/bin/")
SET(BUILD_INPUT_PATH "${EXECUTABLE_OUTPUT_PATH}/input")
SET(BUILD_RELOAD_PATH "${EXECUTABLE_OUTPUT_PATH}/reload")

file(MAKE_DIRECTORY ${BUILD_INPUT_PATH})

add_executable(${PROJECT_NAME})
aux_source_directory(. DIR_SRCS)
target_sources(${PROJECT_NAME} PRIVATE ${DIR_SRCS})
target_link_libraries(${PROJECT_NAME} sphinxsys_3d)
set_target_properties(${PROJECT_NAME} PROPERTIES VS_DEBUGGER_WORKING_DIRECTORY "${EXECUTABLE_OUTPUT_PATH}")

add_test(NAME ${PROJECT_NAME}
    COMMAND ${PROJECT_NAME} --state_recording=${TEST_STATE_RECORDING}
    WORKING_DIRECTORY ${EXECUTABLE_OUTPUT_PATH})

# Fast regressions run in separate directories and cannot overwrite production output.
foreach(case IN ITEMS smoke_liquid smoke_dry smoke_brittle smoke_fracture)
  file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/${case}")
  add_test(NAME peri_capsule_${case} COMMAND ${PROJECT_NAME}
    --case-config "${CMAKE_CURRENT_SOURCE_DIR}/configs/${case}.cfg" --state_recording=false)
  set_tests_properties(peri_capsule_${case} PROPERTIES
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/${case}" TIMEOUT 180)
endforeach()

find_package(Python3 COMPONENTS Interpreter REQUIRED)
foreach(case IN ITEMS smoke_liquid smoke_dry smoke_brittle smoke_fracture)
  set_tests_properties(peri_capsule_${case} PROPERTIES FIXTURES_SETUP ${case}_complete LABELS peri_fast)
  set(damage_check)
  if(case STREQUAL "smoke_fracture")
    set(damage_check --require-damage)
  endif()
  add_test(NAME peri_capsule_${case}_verify COMMAND ${Python3_EXECUTABLE}
    "${CMAKE_SOURCE_DIR}/tools/peri_run.py" verify --run-dir "${CMAKE_CURRENT_BINARY_DIR}/${case}" ${damage_check})
  set_tests_properties(peri_capsule_${case}_verify PROPERTIES FIXTURES_REQUIRED ${case}_complete LABELS peri_fast)
endforeach()

option(PERI_ENABLE_LONG_TESTS "Enable the full-resolution default capsule CTest" OFF)
set_tests_properties(${PROJECT_NAME} PROPERTIES LABELS peri_long)
if(NOT PERI_ENABLE_LONG_TESTS)
  set_tests_properties(${PROJECT_NAME} PROPERTIES DISABLED TRUE)
endif()
add_test(NAME peri_capsule_injected_nan COMMAND ${Python3_EXECUTABLE}
  "${CMAKE_SOURCE_DIR}/tools/peri_run.py" fault-test
  --executable "$<TARGET_FILE:${PROJECT_NAME}>"
  --config "${CMAKE_CURRENT_SOURCE_DIR}/configs/inject_failure.cfg"
  --run-root "${CMAKE_CURRENT_BINARY_DIR}/fault_runs")
set_tests_properties(peri_capsule_injected_nan PROPERTIES LABELS peri_fast TIMEOUT 30)
