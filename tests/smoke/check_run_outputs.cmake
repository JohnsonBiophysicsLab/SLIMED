# Verification step for the executable smoke tests (see tests/CMakeLists.txt).
#
#   cmake -DSMOKE_DIR=<dir> -DSMOKE_FILES=<a.csv,b.csv,...> [-DSMOKE_CHECK_ENERGY=ON] -P check_run_outputs.cmake
#
# Fails unless every listed file exists in SMOKE_DIR and is non-empty and,
# with SMOKE_CHECK_ENERGY, unless EnergyForce.csv has data rows of eight
# finite columns. The multithreading drivers leave SMOKE_CHECK_ENERGY off:
# all their replicas write the one unprefixed EnergyForce.csv at once.

foreach(var SMOKE_DIR SMOKE_FILES)
    if(NOT DEFINED ${var} OR "${${var}}" STREQUAL "")
        message(FATAL_ERROR "check_run_outputs.cmake: ${var} is not set")
    endif()
endforeach()

set(missing "")
string(REPLACE "," ";" files "${SMOKE_FILES}")
foreach(f IN LISTS files)
    set(path "${SMOKE_DIR}/${f}")
    if(NOT EXISTS "${path}")
        list(APPEND missing "${f} (missing)")
    else()
        file(SIZE "${path}" size)
        if(size EQUAL 0)
            list(APPEND missing "${f} (empty)")
        endif()
    endif()
endforeach()
if(missing)
    string(REPLACE ";" ", " missing "${missing}")
    message(FATAL_ERROR "smoke run in ${SMOKE_DIR} did not produce: ${missing}")
endif()

set(energyForce "${SMOKE_DIR}/EnergyForce.csv")
if(SMOKE_CHECK_ENERGY)
    if(NOT EXISTS "${energyForce}")
        message(FATAL_ERROR "smoke run in ${SMOKE_DIR} did not produce EnergyForce.csv")
    endif()
    file(STRINGS "${energyForce}" lines)
    list(LENGTH lines nLines)
    if(nLines LESS 2)
        message(FATAL_ERROR "${energyForce} has no data rows")
    endif()
    list(REMOVE_AT lines 0) # the header
    foreach(line IN LISTS lines)
        string(REPLACE "," ";" cells "${line}")
        list(LENGTH cells nCells)
        if(NOT nCells EQUAL 8)
            message(FATAL_ERROR "${energyForce}: expected 8 columns, got ${nCells}: ${line}")
        endif()
        string(TOLOWER "${line}" lower)
        if(lower MATCHES "nan|inf")
            message(FATAL_ERROR "${energyForce}: non-finite value in row: ${line}")
        endif()
    endforeach()
endif()

message(STATUS "smoke run in ${SMOKE_DIR}: all ${SMOKE_FILES} present")
