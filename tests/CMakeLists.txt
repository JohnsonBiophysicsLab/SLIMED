# GoogleTest unit tests -- the CMake equivalent of `make test`.
#
# Run them with `ctest --test-dir build`, or directly:
#
#     ./build/bin/test_main       # must be run from the repository root
#
# Several tests open data files through paths relative to the working
# directory ("input.params", "./data/example/..."), so ctest is told to run the
# binary from ${PROJECT_SOURCE_DIR} rather than from the build tree.

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
