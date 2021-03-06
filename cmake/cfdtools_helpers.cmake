function(cfdtools_add_library libname)
    # Helper to create properly namespaced library targets
    set(fullname "cfdtools_${libname}")
    add_library(${fullname} ${ARGN})
    add_library("cfdtools::${libname}" ALIAS ${fullname})
    target_link_libraries(${fullname} PRIVATE cfdtools_global)
    target_include_directories(${fullname} PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
    set_target_properties(${fullname} PROPERTIES
        EXPORT_NAME ${libname})
    install(TARGETS ${fullname}
        EXPORT cfdtools
        DESTINATION ${CMAKE_INSTALL_LIBDIR})
    install(FILES README
        DESTINATION "${CMAKE_INSTALL_DATADIR}/cfdtools/doc"
        RENAME ${libname}
        OPTIONAL)
endfunction()

function(cfdtools_add_executable exename)
    # Helper to construct properly namespaced executable targets
    set(fullname "cfdtools_${exename}")
    add_executable(${fullname} ${ARGN})
    install(PROGRAMS $<TARGET_FILE:${fullname}>
        DESTINATION ${CMAKE_INSTALL_BINDIR}
        RENAME ${exename})
    install(FILES README
        DESTINATION "${CMAKE_INSTALL_DATADIR}/cfdtools/doc"
        RENAME "${exename}"
        OPTIONAL)
endfunction()