# set type specific output defaults
include(GNUInstallDirs)

SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}")
SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}")
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}")
SET(CMAKE_Fortran_MODULE_DIRECTORY "${PROJECT_BINARY_DIR}/modules")

set(LIB amos)
add_library(${LIB}
        amos.f90
        d1mach.f90
        dgamln.f
        dsclmr.f
        fdump.f
        i1mach.f90
        r1mach.f90
        zabs.f
        zacai.f
        zacon.f
        zairy.f
        zasyi.f
        zbesh.f
        zbesi.f
        zbesj.f
        zbesk.f
        zbesy.f
        zbinu.f
        zbiry.f
        zbknu.f
        zbuni.f
        zbunk.f
        zdiv.f
        zexp.f
        zkscl.f
        zlog.f
        zmlri.f
        zmlt.f
        zrati.f
        zs1s2.f
        zseri.f
        zshch.f
        zsqrt.f
        zuchk.f
        zunhj.f
        zuni1.f
        zuni2.f
        zunik.f
        zunk1.f
        zunk2.f
        zuoik.f
        zwrsk.f
        )
add_library(${NAMESPACE}${LIB} ALIAS ${LIB})

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
