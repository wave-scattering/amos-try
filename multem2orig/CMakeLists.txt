cmake_minimum_required(VERSION 3.0.2)


SET (PROJECTNAME "multem2orig")

PROJECT (${PROJECTNAME} Fortran)

ENABLE_LANGUAGE(C) 

#FIND_PACKAGE(LAPACK REQUIRED)

#FIND_PACKAGE(HDF5 COMPONENTS Fortran REQUIRED)
#INCLUDE_DIRECTORIES(${HDF5_INCLUDE_DIR})


#SET (CMAKE_Fortran_FLAGS "-O3 -cpp -mtune=native -fopenmp")
#set(CMAKE_Fortran_FLAGS   " -g -fbounds-check -fcheck-array-temporaries -fbacktrace")
set(CMAKE_Fortran_FLAGS   " -g -fcheck-array-temporaries -fbacktrace")

ADD_EXECUTABLE(${PROJECTNAME} "multem_mod.f")



target_link_libraries(${PROJECTNAME} #${BLAS_LIBRARIES}
        #${LAPACK_LIBRARIES} ${HDF5_Fortran_LIBRARIES}
        )

add_custom_command(
        TARGET ${PROJECTNAME} POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
        #        ${CMAKE_SOURCE_DIR}/multem2mod/input2_band.txt
        ${CMAKE_SOURCE_DIR}/multem2mod/input_inoue.txt
#        ${CMAKE_SOURCE_DIR}/multem2mod/fort.10

        ${CMAKE_CURRENT_BINARY_DIR}/fort.10)

