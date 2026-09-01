# CMake build for SLIMED (Subdivision-LImit MEmbrane Dynamics).
#
# This mirrors the hand-written Makefile.  The Makefile selects a build variant
# through the goal name (make serial | omp | dyna | dyna_omp | multi |
# dyna_multi | test); CMake instead configures one tree and builds every
# executable from it, with the parallelisation variant chosen by an option:
#
#     cmake -S . -B build                              # == make serial / dyna
#     cmake -S . -B build -DSLIMED_ENABLE_OPENMP=ON    # == make omp / dyna_omp
#     cmake --build build -j
#     ctest --test-dir build
#
# The multithreading executables (make multi / dyna_multi) are always built;
# they are ordinary std::thread programs and need no special configuration.

cmake_minimum_required(VERSION 3.16)

project(SLIMED
    VERSION 0.1.0
    DESCRIPTION "Subdivision-limit membrane dynamics model"
    HOMEPAGE_URL "https://github.com/JohnsonBiophysicsLab/SLIMED"
    LANGUAGES CXX)

#---------------------------------------------------------------------------
# Options
#---------------------------------------------------------------------------

option(SLIMED_ENABLE_OPENMP   "Build with OpenMP parallelisation (make omp)" OFF)
option(SLIMED_ENABLE_COVERAGE "Instrument the build for gcov/lcov coverage"  OFF)
option(SLIMED_BUILD_TESTS     "Build the GoogleTest unit tests"              ON)
option(SLIMED_ENABLE_CUDA     "Build the CUDA force backend (needs nvcc)"    OFF)
set(SLIMED_CUDA_ARCHITECTURES "70;80;90" CACHE STRING
    "GPU architectures to generate code for when SLIMED_ENABLE_CUDA is ON")

#---------------------------------------------------------------------------
# Global settings
#---------------------------------------------------------------------------

set(CMAKE_CXX_STANDARD 14)          # matches CXXFLAGS = -std=c++14
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# The Makefile compiles at -O3 unconditionally; Release is the closest match.
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY
        STRINGS Debug Release RelWithDebInfo MinSizeRel)
endif()

# Keep the Makefile's layout: executables land in <build>/bin.
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

# Handy for clangd / IDE tooling.
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

#---------------------------------------------------------------------------
# Dependencies
#---------------------------------------------------------------------------

# GSL is used by Parameters.hpp, Linear_algebra.hpp and Scaffolding_points.cpp.
find_package(GSL REQUIRED)

# std::thread in the *_multithreading executables, and GoogleTest, need pthread.
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

# Not REQUIRED: a serial build only needs omp.h to be reachable (see below).
find_package(OpenMP)

if(APPLE AND NOT OpenMP_CXX_FOUND)
    # Homebrew keeps libomp keg-only, so it sits outside Apple Clang's default
    # search path and the probe above misses it.  Retry with the cellar prefix.
    find_program(SLIMED_BREW_EXECUTABLE brew)
    if(SLIMED_BREW_EXECUTABLE)
        execute_process(
            COMMAND "${SLIMED_BREW_EXECUTABLE}" --prefix libomp
            OUTPUT_VARIABLE SLIMED_LIBOMP_PREFIX
            OUTPUT_STRIP_TRAILING_WHITESPACE
            ERROR_QUIET)
        if(EXISTS "${SLIMED_LIBOMP_PREFIX}/include/omp.h")
            set(OpenMP_ROOT "${SLIMED_LIBOMP_PREFIX}")
            find_package(OpenMP)
        endif()
    endif()
endif()

#---------------------------------------------------------------------------
# Coverage
#---------------------------------------------------------------------------

# Carries the gcov instrumentation flags; equivalent to `make ... COVERAGE=1`.
add_library(slimed_coverage INTERFACE)

if(SLIMED_ENABLE_COVERAGE)
    if(NOT CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        message(FATAL_ERROR
            "SLIMED_ENABLE_COVERAGE needs GCC or Clang, not "
            "'${CMAKE_CXX_COMPILER_ID}'.")
    endif()
    target_compile_options(slimed_coverage INTERFACE -O0 -g --coverage)
    target_link_options(slimed_coverage INTERFACE --coverage)
endif()

#---------------------------------------------------------------------------
# Core library
#---------------------------------------------------------------------------

# The Makefile discovers sources with `find $(SDIR) -name "*.cpp"`, so globbing
# keeps the two build systems in step.  CONFIGURE_DEPENDS re-runs the glob when
# the source tree changes.
file(GLOB_RECURSE SLIMED_SOURCES CONFIGURE_DEPENDS "${PROJECT_SOURCE_DIR}/src/*.cpp")

if(NOT SLIMED_SOURCES)
    message(FATAL_ERROR "No sources found under ${PROJECT_SOURCE_DIR}/src.")
endif()

# The CUDA backend and its stub define the same symbols; exactly one of them is
# compiled.  The stub is what a build without a GPU toolchain links, so that
# callers can be written as "try the device, fall back to the host" with no
# #ifdef at the call site.
if(SLIMED_ENABLE_CUDA)
    list(REMOVE_ITEM SLIMED_SOURCES
        "${PROJECT_SOURCE_DIR}/src/cuda/Cuda_force_backend_stub.cpp")
    list(APPEND SLIMED_SOURCES "${PROJECT_SOURCE_DIR}/src/cuda/Cuda_force_backend.cu")
endif()

add_library(slimed_core STATIC ${SLIMED_SOURCES})
add_library(SLIMED::core ALIAS slimed_core)

target_include_directories(slimed_core PUBLIC "${PROJECT_SOURCE_DIR}/include")

target_link_libraries(slimed_core
    PUBLIC  GSL::gsl GSL::gslcblas Threads::Threads
    PRIVATE slimed_coverage)

#---------------------------------------------------------------------------
# OpenMP
#---------------------------------------------------------------------------
#
# -fopenmp and -DOMP must be enabled together, never one alone.  Most `#pragma
# omp parallel for` directives in src/ are unguarded, but the loop bodies that
# need a thread id -- e.g. Compute_energy_and_force_on_mesh.cpp:50 -- read it
# from omp_get_thread_num() only under `#ifdef OMP`, and otherwise hard-code 0.
# Compiling the pragmas without the define would let every thread write to the
# same per-thread accumulation buffer.

if(SLIMED_ENABLE_OPENMP)
    if(NOT OpenMP_CXX_FOUND)
        message(FATAL_ERROR
            "SLIMED_ENABLE_OPENMP=ON but no OpenMP runtime was found. "
            "Install libomp-dev (Linux) or `brew install libomp` (macOS), or "
            "configure with -DSLIMED_ENABLE_OPENMP=OFF.")
    endif()
    target_link_libraries(slimed_core PUBLIC OpenMP::OpenMP_CXX)
    target_compile_definitions(slimed_core PUBLIC OMP)
elseif(OpenMP_CXX_FOUND)
    # Serial build.  Mesh.hpp, Dynamics.hpp, Force.hpp and Gauss_quadrature.hpp
    # `#include <omp.h>` unconditionally, so the header search path is still
    # required -- but deliberately not the -fopenmp flag (see above).
    target_include_directories(slimed_core SYSTEM PUBLIC ${OpenMP_CXX_INCLUDE_DIRS})
else()
    message(WARNING
        "OpenMP was not found.  Several headers include <omp.h> unconditionally, "
        "so the build will fail unless omp.h is already on the compiler's "
        "include path (it is with GCC; it is not with Apple Clang).")
endif()

#---------------------------------------------------------------------------
# CUDA
#---------------------------------------------------------------------------
#
# Off by default: the CPU path is complete on its own, and most machines that
# build this have no NVIDIA GPU.  When it is on, the only extra source is
# src/cuda/Cuda_force_backend.cu -- the kernel bodies are the header-only,
# __host__ __device__ functions in include/energy_force/Patch_kernel.inl and
# include/cuda/Force_kernels.inl that the CPU path already runs and the tests
# already pin against the production force evaluation.
#
# Verify a new GPU build before trusting it:
#
#     ctest --test-dir build -R CudaForceBackend
#
# CudaForceBackendTest.MatchesTheProductionForceEvaluation runs both paths over
# the same meshes and compares every energy and force.  It skips, rather than
# passes, when no device is visible.

if(SLIMED_ENABLE_CUDA)
    include(CheckLanguage)
    check_language(CUDA)
    if(NOT CMAKE_CUDA_COMPILER)
        message(FATAL_ERROR
            "SLIMED_ENABLE_CUDA=ON but no CUDA compiler was found.  Put nvcc on PATH "
            "(module load cuda, on most clusters), point CMAKE_CUDA_COMPILER at it, "
            "or configure with -DSLIMED_ENABLE_CUDA=OFF.")
    endif()
    enable_language(CUDA)

    set_target_properties(slimed_core PROPERTIES
        CUDA_ARCHITECTURES "${SLIMED_CUDA_ARCHITECTURES}"
        # The kernel bodies are inline and defined in headers, so each
        # translation unit carries its own device code and no cross-TU device
        # linking is needed.
        CUDA_SEPARABLE_COMPILATION OFF)

    # Mesh.hpp and the headers it pulls in are ordinary host C++; nvcc passes
    # them to the host compiler, which needs the same OpenMP include path the
    # C++ sources get.
    if(OpenMP_CXX_FOUND)
        target_include_directories(slimed_core SYSTEM PUBLIC ${OpenMP_CXX_INCLUDE_DIRS})
    endif()

    target_compile_definitions(slimed_core PUBLIC SLIMED_WITH_CUDA)
    message(STATUS "CUDA force backend enabled for architectures: ${SLIMED_CUDA_ARCHITECTURES}")
endif()

#---------------------------------------------------------------------------
# Executables
#---------------------------------------------------------------------------

# Each name is both the target and its EXEs/<name>.cpp entry point.
set(SLIMED_APPS
    continuum_membrane                  # make serial / make omp
    continuum_membrane_multithreading   # make multi
    membrane_dynamics                   # make dyna / make dyna_omp
    membrane_dynamics_multithreading)   # make dyna_multi

foreach(app IN LISTS SLIMED_APPS)
    add_executable(${app} "${PROJECT_SOURCE_DIR}/EXEs/${app}.cpp")
    target_link_libraries(${app} PRIVATE slimed_core slimed_coverage)
endforeach()

#---------------------------------------------------------------------------
# Tests
#---------------------------------------------------------------------------

enable_testing()

if(SLIMED_BUILD_TESTS)
    add_subdirectory(tests)
endif()

#---------------------------------------------------------------------------
# Configuration summary
#---------------------------------------------------------------------------

message(STATUS "")
message(STATUS "  ${PROJECT_NAME} ${PROJECT_VERSION}")
message(STATUS "    C++ compiler ...... ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "    Build type ........ ${CMAKE_BUILD_TYPE}")
message(STATUS "    GSL ............... ${GSL_VERSION}")
message(STATUS "    OpenMP (-DOMP) .... ${SLIMED_ENABLE_OPENMP}")
message(STATUS "    Coverage .......... ${SLIMED_ENABLE_COVERAGE}")
message(STATUS "    Unit tests ........ ${SLIMED_BUILD_TESTS}")
message(STATUS "")
