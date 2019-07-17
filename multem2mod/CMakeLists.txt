#cmake_minimum_required(VERSION 3.0.2)


SET (PROJECTNAME "multem2")

#PROJECT (${PROJECTNAME} Fortran)

#ENABLE_LANGUAGE(C)

#FIND_PACKAGE(LAPACK REQUIRED)

#FIND_PACKAGE(HDF5 COMPONENTS Fortran REQUIRED)
#INCLUDE_DIRECTORIES(${HDF5_INCLUDE_DIR})


#SET (CMAKE_Fortran_FLAGS "-O3 -cpp -mtune=native -fopenmp")

ADD_EXECUTABLE(${PROJECTNAME} "multem_mod.f"
        "bessel_switch.f"
        )



target_link_libraries(${PROJECTNAME} amos ${BLAS_LIBRARIES}
        ${LAPACK_LIBRARIES} #${HDF5_Fortran_LIBRARIES}
        )


