cfdtools_add_executable(thin_grid_2d
  thin_grid_2d.f90
)
target_link_libraries(cfdtools_thin_grid_2d
  cfdtools_xyzq_io
)
