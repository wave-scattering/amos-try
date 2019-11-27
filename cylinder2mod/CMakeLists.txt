#cmake_minimum_required(VERSION 3.0.2)
SET (PROJECTNAME "cylinder2mod")
PROJECT (${PROJECTNAME} Fortran)
FIND_PACKAGE(LAPACK REQUIRED)

include_directories(${CMAKE_BINARY_DIR}/amos
        ${CMAKE_BINARY_DIR}/lapack_wrapper
        ${CMAKE_BINARY_DIR}/errfun
        )
#include_directories(${CMAKE_BINARY_DIR}/errfun
#        )

add_library(constants constants.f90)

add_library(model_parameters model_parameters.f90)
target_link_libraries(model_parameters constants FiNeR::FiNeR)

add_library(special_functions special_functions.f90)
target_link_libraries(special_functions constants)

add_library(cylinder_blas cylinder_blas.f90)
target_link_libraries(cylinder_blas constants)

add_library(libcylinder libcylinder.f90)
target_link_libraries(libcylinder
        constants
        special_functions
        model_parameters
        cylinder_blas
        )

ADD_EXECUTABLE(${PROJECTNAME}
        tmtaxsp1.f
        "axspartclnew4.f"
        gnzbess.f
        biga.f
        chebev.f # TODO replace with non NR
        gnricbessh.f
        beschb.f # TODO replace with non NR
        sordalc.f
#        lpd.f
        )
target_link_libraries(${PROJECTNAME}
        constants
        model_parameters
        libcylinder
        errfun
        amos
        dense_solve
        ${BLAS_LIBRARIES}
        ${LAPACK_LIBRARIES}
        )

#add_custom_command(
#        TARGET ${PROJECTNAME} POST_BUILD
#        COMMAND ${CMAKE_COMMAND} -E copy
#        ${CMAKE_SOURCE_DIR}/cylinder2mod/default.ini
#        ${CMAKE_CURRENT_BINARY_DIR}/fort.10)



