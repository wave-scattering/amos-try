#cmake_minimum_required(VERSION 3.0.2)
SET (PROJECTNAME "mishchenko")
PROJECT (${PROJECTNAME} Fortran)

#include_directories(${CMAKE_BINARY_DIR}/amos
#        ${CMAKE_BINARY_DIR}/lapack_wrapper
#        ${CMAKE_BINARY_DIR}/errfun
#        )
#include_directories(${CMAKE_BINARY_DIR}/errfun
#        )
ADD_EXECUTABLE(${PROJECTNAME} ampld.lp.f
        lpd.f
        )
#target_link_libraries(${PROJECTNAME}
#        errfun
#        amos
#        dense_solve
#        )


