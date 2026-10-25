# GoogleTest unit tests -- the CMake equivalent of `make test`.
#
# Run them with `ctest --test-dir build`, or directly:
#
#     ./build/bin/test_main       # must be run from the repository root
#
# Several tests open data files through paths relative to the working
# directory ("input.params", "./data/example/..."), so ctest is told to run the
# binary from ${PROJECT_SOURCE_DIR} rather than from the build tree.

#---------------------------------------------------------------------------
# Executable smoke tests
#---------------------------------------------------------------------------
#
# The four simulation executables take no arguments: they read fixed-name
# params files from their working directory (input.params; input1/input2.params
# for continuum_membrane_multithreading; input1 ... input8.params for
# membrane_dynamics_multithreading) and write fixed-name outputs into it. Each smoke test runs
# one of them in its own scratch directory under the build tree, on a tiny
# sheet from tests/smoke/, in three ctest steps chained by fixtures:
#
#   smoke_<app>_setup   empty the directory, copy the params template in
#   smoke_<app>         run the executable there
#   smoke_<app>_check   the outputs it promises exist (and are finite)
#
# A missing input file would silently fall back to the full-size defaults, so
# the multithreading drivers get every one of their input files. These do not
# need GoogleTest and are defined before the GTest check below.

set(SLIMED_SMOKE_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/smoke")
set(SLIMED_SMOKE_RUN_ROOT "${CMAKE_CURRENT_BINARY_DIR}/smoke_runs")

# slimed_add_smoke_test(<app> <template> <comma-separated input names>
#                       <comma-separated expected outputs> <check EnergyForce.csv: ON|OFF>)
function(slimed_add_smoke_test app template names outputs checkEnergy)
    set(runDir "${SLIMED_SMOKE_RUN_ROOT}/${app}")

    add_test(NAME smoke_${app}_setup
        COMMAND "${CMAKE_COMMAND}"
            "-DSMOKE_DIR=${runDir}"
            "-DSMOKE_TEMPLATE=${SLIMED_SMOKE_SOURCE_DIR}/${template}"
            "-DSMOKE_NAMES=${names}"
            -P "${SLIMED_SMOKE_SOURCE_DIR}/prepare_run_dir.cmake")
    set_tests_properties(smoke_${app}_setup PROPERTIES
        FIXTURES_SETUP smoke_${app}_dir
        LABELS smoke
        TIMEOUT 30)

    add_test(NAME smoke_${app} COMMAND $<TARGET_FILE:${app}>)
    set_tests_properties(smoke_${app} PROPERTIES
        WORKING_DIRECTORY "${runDir}"
        FIXTURES_REQUIRED smoke_${app}_dir
        FIXTURES_SETUP smoke_${app}_ran
        LABELS smoke
        TIMEOUT 60)

    add_test(NAME smoke_${app}_check
        COMMAND "${CMAKE_COMMAND}"
            "-DSMOKE_DIR=${runDir}"
            "-DSMOKE_FILES=${outputs}"
            "-DSMOKE_CHECK_ENERGY=${checkEnergy}"
            -P "${SLIMED_SMOKE_SOURCE_DIR}/check_run_outputs.cmake")
    set_tests_properties(smoke_${app}_check PROPERTIES
        FIXTURES_REQUIRED "smoke_${app}_dir;smoke_${app}_ran"
        LABELS smoke
        TIMEOUT 30)
endfunction()

slimed_add_smoke_test(continuum_membrane
    continuum_membrane.params
    "input"
    "face.csv,vertex_begin.csv,vertex_initial_step_0.csv,EnergyForce.csv,ElementFaceEnergy.csv,vertexfinal.csv,smoke_restart.chk"
    ON)

# Two replicas, input1 and input2; run_flat's output names are not prefixed,
# so both write the same files.
slimed_add_smoke_test(continuum_membrane_multithreading
    continuum_membrane.params
    "input1,input2"
    "face.csv,vertex_begin.csv,EnergyForce.csv,ElementFaceEnergy.csv,vertexfinal.csv"
    OFF)

slimed_add_smoke_test(membrane_dynamics
    membrane_dynamics_trajectory.params
    "input"
    "inputface.csv,inputvertex_begin.csv,inputvertex_type_begin.csv,inputvertex_final.csv,inputvertex_type_final.csv,meshpointinput.csv,surfacepointinput.csv,faceenergyinput.csv,vertexforceinput.csv,EnergyForce.csv"
    ON)

set(SLIMED_SMOKE_DYNAMICS_REPLICAS "")
set(SLIMED_SMOKE_DYNAMICS_REPLICA_OUTPUTS "")
foreach(i RANGE 1 8)
    list(APPEND SLIMED_SMOKE_DYNAMICS_REPLICAS "input${i}")
    list(APPEND SLIMED_SMOKE_DYNAMICS_REPLICA_OUTPUTS
        "input${i}face.csv" "input${i}vertex_begin.csv" "input${i}vertex_final.csv")
endforeach()
list(APPEND SLIMED_SMOKE_DYNAMICS_REPLICA_OUTPUTS "EnergyForce.csv")
string(REPLACE ";" "," SLIMED_SMOKE_DYNAMICS_REPLICAS "${SLIMED_SMOKE_DYNAMICS_REPLICAS}")
string(REPLACE ";" "," SLIMED_SMOKE_DYNAMICS_REPLICA_OUTPUTS "${SLIMED_SMOKE_DYNAMICS_REPLICA_OUTPUTS}")

slimed_add_smoke_test(membrane_dynamics_multithreading
    membrane_dynamics.params
    "${SLIMED_SMOKE_DYNAMICS_REPLICAS}"
    "${SLIMED_SMOKE_DYNAMICS_REPLICA_OUTPUTS}"
    OFF)

#---------------------------------------------------------------------------
# GoogleTest unit tests
#---------------------------------------------------------------------------

find_package(GTest QUIET)

if(NOT GTest_FOUND)
    message(WARNING
        "GoogleTest was not found -- skipping the unit tests.  Install it with "
        "`sudo apt-get install libgtest-dev` (Linux) or `brew install googletest` "
        "(macOS), or configure with -DSLIMED_BUILD_TESTS=OFF to silence this.")
    return()
endif()

file(GLOB SLIMED_TEST_SOURCES CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")

# EXEs/test_main.cpp supplies its own main(), so link gtest rather than
# gtest_main.
add_executable(test_main
    "${PROJECT_SOURCE_DIR}/EXEs/test_main.cpp"
    ${SLIMED_TEST_SOURCES})

target_include_directories(test_main PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")

target_link_libraries(test_main
    PRIVATE slimed_core GTest::gtest Threads::Threads slimed_coverage)

include(GoogleTest)
gtest_discover_tests(test_main
    WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
    PROPERTIES WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}")
