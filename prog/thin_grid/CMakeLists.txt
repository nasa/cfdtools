cfdtools_add_executable(thin_grid
  thin_grid.f90
)
target_link_libraries(cfdtools_thin_grid
  cfdtools_xyzq_io
)
