#cmake_minimum_required(VERSION 3.0.2)


SET (PROJECTNAME "quadpack_test")
PROJECT (${PROJECTNAME} Fortran)

ADD_EXECUTABLE(${PROJECTNAME}
        quadpack_double_prb.f90
        )
target_link_libraries(${PROJECTNAME}
        quadpack::quadpack
        )


