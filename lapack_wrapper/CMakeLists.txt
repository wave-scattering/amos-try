project(dense_solve)
set(SRC
    ${CMAKE_CURRENT_SOURCE_DIR}/dense_solve.f90
)

ADD_LIBRARY(dense_solve ${SRC})
