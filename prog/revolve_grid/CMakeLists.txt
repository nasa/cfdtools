cfdtools_add_executable(revolve_grid
  revolve_grid.f90
)
target_link_libraries(cfdtools_revolve_grid
  cfdtools_geomlib
  cfdtools_gridlib
  cfdtools_xyzq_io
)
