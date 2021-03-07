cfdtools_add_library(xyzq_io
  f_io.f90
  xyq_io.f90
  xyzq_io.f90
)
target_link_libraries(cfdtools_xyzq_io PUBLIC cfdtools_grid_block_structure)
