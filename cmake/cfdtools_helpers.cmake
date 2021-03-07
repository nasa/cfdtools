#------------------------------------------------------------------------------
# Boilerplate for library targets
#------------------------------------------------------------------------------
function(cfdtools_add_library libname)
    set(fullname "cfdtools_${libname}")
    set(doc_dir "${CMAKE_INSTALL_DATADIR}/cfdtools/doc")
    set(mod_build_dir "${CMAKE_CURRENT_BINARY_DIR}/modules")
    set(mod_install_dir "${CMAKE_INSTALL_INCLUDEDIR}/cfdtools/${libname}")

    # Basic definition
    add_library(${fullname} ${ARGN})
    add_library("cfdtools::${libname}" ALIAS ${fullname})
    set_target_properties(${fullname} PROPERTIES
        EXPORT_NAME ${libname}
        Fortran_MODULE_DIRECTORY ${mod_build_dir}
    )

    # Default dependencies
    target_link_libraries(${fullname} PRIVATE cfdtools_global)

    # Expose module files
    target_include_directories(${fullname} PUBLIC
        $<BUILD_INTERFACE:${mod_build_dir}>
        $<INSTALL_INTERFACE:${mod_install_dir}>
    )

    # Installation
    install(TARGETS ${fullname}
        EXPORT cfdtools
        DESTINATION ${CMAKE_INSTALL_LIBDIR})
    install(DIRECTORY "${mod_build_dir}/"  # Trailing slash => copy contents
        DESTINATION ${mod_install_dir})
    install(FILES README
        DESTINATION ${doc_dir}
        RENAME ${libname}
        OPTIONAL)

endfunction()

#------------------------------------------------------------------------------
# Boilerplate for executable targets
#------------------------------------------------------------------------------
function(cfdtools_add_executable exename)
    set(fullname "cfdtools_${exename}")
    add_executable(${fullname} ${ARGN})
    target_link_libraries(${fullname} cfdtools_global)
    install(PROGRAMS $<TARGET_FILE:${fullname}>
        DESTINATION ${CMAKE_INSTALL_BINDIR}
        RENAME ${exename})
    install(FILES README
        DESTINATION "${CMAKE_INSTALL_DATADIR}/cfdtools/doc"
        RENAME "${exename}"
        OPTIONAL)
endfunction()