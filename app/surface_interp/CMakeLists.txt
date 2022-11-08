cfdtools_add_executable(surface_interp surface_interp.f90)
target_link_libraries(cfdtools_surface_interp
    cfdtools_tecplot_adt
    cfdtools_eigenlib
    cfdtools_geomlib
    cfdtools_linsys
    cfdtools_progtools
    cfdtools_tecplot_gridlib
    cfdtools_tecplot_io
    cfdtools_tecplot_xyzq_io
    cfdtools_triangulation_io
    cfdtools_ugridlib
)
