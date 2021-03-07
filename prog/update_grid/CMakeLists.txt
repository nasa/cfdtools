cfdtools_add_executable(update_grid
  update_grid.f90
)
target_link_libraries(cfdtools_update_grid
  cfdtools_progtools
  cfdtools_geomlib
  cfdtools_xyzq_io
)
