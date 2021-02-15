# set type specific output defaults
include(GNUInstallDirs)


SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}")
SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}")
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}")
SET(CMAKE_Fortran_MODULE_DIRECTORY "${PROJECT_BINARY_DIR}/modules")

set(LIB dense_solve)
add_library(${LIB}
        dense_solve.f90
        dense_solve_impl.f90
        )
add_library(${NAMESPACE}${LIB} ALIAS ${LIB})
target_link_libraries(dense_solve ${BLAS_LIBRARIES}
        ${LAPACK_LIBRARIES} #${HDF5_Fortran_LIBRARIES}
        )
target_include_directories(${LIB}
        INTERFACE
        $<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>
        )


#target_compile_definitions(${LIB}
#    PRIVATE
#        ${r16p_supported}
#        ${ascii_supported}
#        ${ascii_neq_default}
#        ${ucs4_supported}
#)

set_target_properties(${LIB} PROPERTIES
        VERSION
        ${PROJECT_VERSION}
        SOVERSION
        ${PROJECT_VERSION_MAJOR}
        )

# installation and export of targets
install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}
        COMPONENT Developement
        )

install(TARGETS ${LIB} EXPORT ${TARGETS_EXPORT_NAME}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT RuntimeLibraries
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT RuntimeLibraries
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}
        )
