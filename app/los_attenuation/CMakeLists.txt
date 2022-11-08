cfdtools_add_executable(los_attenuation los_attenuation.f90)
target_link_libraries(cfdtools_los_attenuation
    cfdtools_tecplot_io
    cfdtools_xyzq_io
    cfdtools_geomlib
    cfdtools_integratelib
    cfdtools_interplib
    cfdtools_searchlib
    cfdtools_numodules
    cfdtools_prmodules
)
