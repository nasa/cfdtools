cfdtools_add_executable(tri_diffs tri_diffs.f90)
target_link_libraries(cfdtools_tri_diffs
    cfdtools_adt
    cfdtools_eigenlib
    cfdtools_geomlib
    cfdtools_linsys
    cfdtools_progtools
    cfdtools_triangulation_io
    cfdtools_ugridlib
    cfdtools_xyzq_io
)
