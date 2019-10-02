project(try-besj Fortran)

ADD_EXECUTABLE(besj-f77 test-besj.f)
target_link_libraries(besj-f77 amos)


include_directories(${CMAKE_BINARY_DIR}/amos)
ADD_EXECUTABLE(besj-f90 test_besselj.f90)
target_link_libraries(besj-f90 amos)


