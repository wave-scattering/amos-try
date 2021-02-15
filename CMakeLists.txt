cmake_minimum_required(VERSION 3.0.2)

# apt install libblas-dev liblapack-dev
SET (PROJECTNAME "top-besj")
PROJECT (${PROJECTNAME} Fortran)
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)


# Make sure that CMAKE_BUILD_TYPE is either Debug or Release:
if (NOT CMAKE_BUILD_TYPE)

#    set(CMAKE_BUILD_TYPE Release
#            CACHE STRING "Build type (Debug, Release)" FORCE)

    set(CMAKE_BUILD_TYPE Debug
            CACHE STRING "Build type (Debug, Release)" FORCE)

endif ()
if (NOT (CMAKE_BUILD_TYPE STREQUAL "Debug" OR
        CMAKE_BUILD_TYPE STREQUAL "Release"))
    message("${CMAKE_BUILD_TYPE}")
    message(FATAL_ERROR "CMAKE_BUILD_TYPE must be one of: Debug, Release (current value: '${CMAKE_BUILD_TYPE}')")
endif ()

if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    execute_process(COMMAND ${CMAKE_Fortran_COMPILER} -dumpversion
            OUTPUT_VARIABLE GCC_VERSION)
    if (GCC_VERSION VERSION_LESS 4.5)
        set(common "-std=f95 -Wall -Werror -Wextra -Wimplicit-interface -fPIC -fmax-errors=1")
        set(CMAKE_Fortran_FLAGS_RELEASE "${common} -O3 -march=native -ffast-math -funroll-loops")
        set(CMAKE_Fortran_FLAGS_DEBUG   "${common} -g -fbounds-check -fcheck-array-temporaries -fbacktrace")
    endif()

    # gfortran
    # Enable this if you want to check for single/double corruption (and use
    # the Debug build):
#    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -fdefault-real-8")
endif ()

#set(CMAKE_BUILD_TYPE Release
#        CACHE STRING "Build type (Debug, Release)" FORCE)
#set(CMAKE_Fortran_FLAGS_RELEASE "-std=f2003 -O3 -march=native -ffast-math -funroll-loops -mtune=native -fopenmp" )
add_subdirectory(amos)
#add_subdirectory(try-zbesj)
add_subdirectory(lapack_wrapper)
#add_subdirectory(erf)
add_subdirectory(erf_MIT)
add_subdirectory(FiNeR)
#add_subdirectory(multem2orig)
#add_subdirectory(cylinder2mod)
add_subdirectory(mishchenko_dp_fixed_particles)
add_subdirectory(quadpack)
set(common "-std=f2008")

#set(CMAKE_BUILD_TYPE Debug
#        CACHE STRING "Build type (Debug, Release)" FORCE)
#set(warn "-Wall -Wunused-variable -Wimplicit-interface -Wunused-dummy-argument -Wextra -Wargument-mismatch  -fPIC -fmax-errors=1")
#set(CMAKE_Fortran_FLAGS_DEBUG   "${common} ${warn} -g -fbounds-check -fcheck-array-temporaries -fbacktrace")

set(CMAKE_BUILD_TYPE Release
            CACHE STRING "Build type (Debug, Release)" FORCE)
#set(CMAKE_Fortran_FLAGS_RELEASE "${common} -O3 -march=native -ffast-math -funroll-loops -mtune=native -fopenmp" )
set(CMAKE_Fortran_FLAGS_RELEASE "${common} -O3 -march=native -ffast-math -funroll-loops -mtune=native" )

#add_subdirectory(cylinder)


add_subdirectory(multem2mod)
add_subdirectory(read-nk)

add_subdirectory(cylinder2mod)
