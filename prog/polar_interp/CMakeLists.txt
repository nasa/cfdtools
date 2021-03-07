cfdtools_add_executable(polar_interp
  polar_interp.f90
)
target_link_libraries(cfdtools_polar_interp
  cfdtools_geomlib
  cfdtools_gridlib
  cfdtools_interplib
  cfdtools_numodules
  cfdtools_searchlib
)

cfdtools_add_executable(reorder_spokes reorder_spokes.f90)
target_link_libraries(cfdtools_reorder_spokes
  cfdtools_table_io
)
