cfdtools_add_executable(neqair_integration
    neqair_integration.f90
)
target_link_libraries(cfdtools_neqair_integration
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

cfdtools_add_executable(cone_angles cone_angles.f90)

cfdtools_add_executable(response_function response_function.f90)
target_link_libraries(cfdtools_response_function cfdtools_interplib cfdtools_table_io)

