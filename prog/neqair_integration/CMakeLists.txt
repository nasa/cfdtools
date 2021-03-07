add_executable(neqair_integration
    neqair_integration.f90
)
target_link_libraries(neqair_integration PRIVATE
    cfdtools_global
    cfdtools_eigenlib
    cfdtools_geomlib
    cfdtools_gridlib
    cfdtools_integratelib
    cfdtools_interplib
    cfdtools_intrinsics
    cfdtools_numodules
    cfdtools_progtools
    cfdtools_searchlib
    cfdtools_triangulation_io
    cfdtools_ugridlib
    cfdtools_xyzq_io
)

add_executable(cone_angles cone_angles.f90)
target_link_libraries(cone_angles PRIVATE cfdtools_global)

add_executable(response_function response_function.f90)
target_link_libraries(response_function PRIVATE cfdtools_global cfdtools_interplib cfdtools_table_io)

install(TARGETS neqair_integration cone_angles response_function
        DESTINATION ${CMAKE_INSTALL_BINDIR})
