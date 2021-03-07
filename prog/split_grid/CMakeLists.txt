cfdtools_add_executable(split_grid
  split_grid.f
)
target_link_libraries(cfdtools_split_grid
  cfdtools_xyzq_io
)
