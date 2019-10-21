#cmake_minimum_required(VERSION 3.0.2)
SET (PROJECTNAME "cylinder2mod")
PROJECT (${PROJECTNAME} Fortran)

include_directories(${CMAKE_BINARY_DIR}/amos
        ${CMAKE_BINARY_DIR}/lapack_wrapper
        ${CMAKE_BINARY_DIR}/errfun
        )
#include_directories(${CMAKE_BINARY_DIR}/errfun
#        )

add_library(libcylinder libcylinder.f90)

ADD_EXECUTABLE(${PROJECTNAME}
        tmtaxsp1.f
        "axspartclnew4.f"
        gnzbess.f
        #bess.f # TODO use amos lib
        bessik.f # TODO use amos lib
        bessjy.f # TODO use amos lib
        biga.f
        chebev.f # TODO replace with non NR
        gnricbessh.f
        beschb.f # TODO replace with non NR
        sordalc.f
        )
target_link_libraries(${PROJECTNAME}
        libcylinder
        errfun
        amos
        dense_solve
        )


