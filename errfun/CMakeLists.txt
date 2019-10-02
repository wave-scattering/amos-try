project(errfun Fortran)
set(SRC
    ${CMAKE_CURRENT_SOURCE_DIR}/errfun.f90
)

ADD_LIBRARY(errfun ${SRC})
