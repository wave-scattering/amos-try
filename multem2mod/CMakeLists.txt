#cmake_minimum_required(VERSION 3.0.2)


SET (PROJECTNAME "multem2")

PROJECT (${PROJECTNAME} Fortran)

#ENABLE_LANGUAGE(C)

#FIND_PACKAGE(LAPACK REQUIRED)

#FIND_PACKAGE(HDF5 COMPONENTS Fortran REQUIRED)
#INCLUDE_DIRECTORIES(${HDF5_INCLUDE_DIR})


#SET (CMAKE_Fortran_FLAGS "-O3 -cpp -mtune=native -fopenmp")

include_directories(${CMAKE_BINARY_DIR}/amos
        ${CMAKE_BINARY_DIR}/lapack_wrapper
        ${CMAKE_BINARY_DIR}/errfun
        )
add_library(multem_blas multem_blas.f90)
target_link_libraries(multem_blas dense_solve)

add_library(libmultem2b libmultem2b.f90)
target_link_libraries(libmultem2b amos dense_solve multem_blas errfun)


ADD_EXECUTABLE(${PROJECTNAME} "multem_mod.f90"
        )
target_link_libraries(${PROJECTNAME} ${BLAS_LIBRARIES}
        dense_solve libmultem2b multem_blas
        errfun
        ${LAPACK_LIBRARIES} #${HDF5_Fortran_LIBRARIES}
        )
add_custom_command(
        TARGET ${PROJECTNAME} POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
        ${CMAKE_SOURCE_DIR}/multem2mod/input2_band.txt
        ${CMAKE_CURRENT_BINARY_DIR}/fort.10)


