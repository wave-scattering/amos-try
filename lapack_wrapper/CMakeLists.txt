project(dense_solve Fortran)
set(SRC
    ${CMAKE_CURRENT_SOURCE_DIR}/dense_solve.f90
)

ADD_LIBRARY(dense_solve ${SRC})
target_link_libraries(dense_solve ${BLAS_LIBRARIES}
        ${LAPACK_LIBRARIES} #${HDF5_Fortran_LIBRARIES}
        )
