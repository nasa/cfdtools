cfdtools_add_executable(uslos uslos.f90)
target_link_libraries(cfdtools_uslos
    cfdtools_adt
    cfdtools_eigenlib
    cfdtools_geomlib
    cfdtools_gridlib
    cfdtools_interp3d
    cfdtools_intrinsics
    cfdtools_linsys
    cfdtools_numodules
    cfdtools_optlib
    cfdtools_progtools
    cfdtools_triangulation_io
    cfdtools_ugridlib
    cfdtools_xyzq_io
)
