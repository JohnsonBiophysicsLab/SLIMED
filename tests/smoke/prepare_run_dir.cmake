# Fixture setup for the executable smoke tests (see tests/CMakeLists.txt).
#
#   cmake -DSMOKE_DIR=<dir> -DSMOKE_TEMPLATE=<file.params> -DSMOKE_NAMES=<a,b,...> -P prepare_run_dir.cmake
#
# Empties SMOKE_DIR and copies SMOKE_TEMPLATE into it once per name, as
# <name>.params: the executables take no arguments and read fixed file names
# from their working directory.

foreach(var SMOKE_DIR SMOKE_TEMPLATE SMOKE_NAMES)
    if(NOT DEFINED ${var} OR "${${var}}" STREQUAL "")
        message(FATAL_ERROR "prepare_run_dir.cmake: ${var} is not set")
    endif()
endforeach()

if(NOT EXISTS "${SMOKE_TEMPLATE}")
    message(FATAL_ERROR "prepare_run_dir.cmake: no template ${SMOKE_TEMPLATE}")
endif()

file(REMOVE_RECURSE "${SMOKE_DIR}")
file(MAKE_DIRECTORY "${SMOKE_DIR}")

string(REPLACE "," ";" names "${SMOKE_NAMES}")
foreach(name IN LISTS names)
    configure_file("${SMOKE_TEMPLATE}" "${SMOKE_DIR}/${name}.params" COPYONLY)
endforeach()
