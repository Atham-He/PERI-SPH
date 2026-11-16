add_executable(peri_sph_liquid_capsule_240um liquid_capsule.cpp)
target_link_libraries(peri_sph_liquid_capsule_240um PRIVATE sphinxsys_3d)
if(BUILD_TESTING)
    add_test(NAME liquid_capsule_ideal_j2
        COMMAND peri_sph_liquid_capsule_240um --material-audit
        "${CMAKE_CURRENT_BINARY_DIR}/ideal_j2_audit.json")
endif()
