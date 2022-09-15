cfdtools_add_executable(usflowinterp usflowinterp.f90)
target_link_libraries(cfdtools_usflowinterp
    cfdtools_xyzq_io
    cfdtools_triangulation_io
    cfdtools_adt
    cfdtools_gridlib
    cfdtools_ugridlib
    cfdtools_geomlib
    cfdtools_linsys
    cfdtools_progtools
    cfdtools_eigenlib
)
