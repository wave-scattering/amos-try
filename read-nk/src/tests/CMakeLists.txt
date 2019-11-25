#cmake_minimum_required(VERSION 3.0.2)


SET (PROJECTNAME "readnk_test")
PROJECT (${PROJECTNAME} Fortran)

ADD_EXECUTABLE(${PROJECTNAME}
        readnk_test.f90
        )
target_link_libraries(${PROJECTNAME}
        readnk::readnk
        )


