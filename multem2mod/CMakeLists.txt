#cmake_minimum_required(VERSION 3.0.2)


SET (PROJECTNAME "multem2")

#PROJECT (${PROJECTNAME} Fortran)

#ENABLE_LANGUAGE(C)

#FIND_PACKAGE(LAPACK REQUIRED)

#FIND_PACKAGE(HDF5 COMPONENTS Fortran REQUIRED)
#INCLUDE_DIRECTORIES(${HDF5_INCLUDE_DIR})


#SET (CMAKE_Fortran_FLAGS "-O3 -cpp -mtune=native -fopenmp")

include_directories(${CMAKE_BINARY_DIR}/amos
        ${CMAKE_BINARY_DIR}/lapack_wrapper
        ${CMAKE_BINARY_DIR}/errfun
        )
add_library(libmultem2b libmultem2b.f90)
target_link_libraries(libmultem2b amos)

ADD_EXECUTABLE(${PROJECTNAME} "multem_mod.f"
        )
target_link_libraries(${PROJECTNAME} ${BLAS_LIBRARIES}
        dense_solve libmultem2b
        errfun
        ${LAPACK_LIBRARIES} #${HDF5_Fortran_LIBRARIES}
        )


