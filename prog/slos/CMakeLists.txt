add_executable(slos slos.f90)
target_link_libraries(slos
  PRIVATE
    cfdtools_global
    cfdtools_adt
    cfdtools_eigenlib
    cfdtools_geomlib
    cfdtools_gridlib
    cfdtools_interp3d
    cfdtools_intrinsics
    cfdtools_linsys
    cfdtools_numodules
    cfdtools_progtools
    cfdtools_optlib
    cfdtools_searchlib
    cfdtools_triangulation_io
    cfdtools_ugridlib
    cfdtools_xyzq_io
)
install(TARGETS slos DESTINATION ${CMAKE_INSTALL_BINDIR})

