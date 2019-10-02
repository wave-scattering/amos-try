#cmake_minimum_required(VERSION 3.0.2)


SET (PROJECTNAME "errfun_test")
PROJECT (${PROJECTNAME} Fortran)


include_directories(${CMAKE_BINARY_DIR}/errfun
        )
ADD_EXECUTABLE(${PROJECTNAME} "errfun_test.f90"
        )
target_link_libraries(${PROJECTNAME}
        errfun
        )


