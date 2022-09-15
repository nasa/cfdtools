cfdtools_add_executable(surface_interp surface_interp.f90)
target_link_libraries(cfdtools_surface_interp
    cfdtools_adt
    cfdtools_eigenlib
    cfdtools_geomlib
    cfdtools_gridlib
    cfdtools_linsys
    cfdtools_progtools
    cfdtools_tecplot_io
    cfdtools_triangulation_io
    cfdtools_ugridlib
    cfdtools_xyzq_io
)
