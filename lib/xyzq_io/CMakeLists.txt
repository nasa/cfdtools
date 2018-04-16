add_library(xyzq_io
  f_io.f90
  grid_block_structure.f90
  xyq_io.f90
  xyq_io.f90.before-xyq_2d_to_3d
  xyzq_io.f90
  xyzq_io.f90.misguided-change
)
target_link_libraries(xyzq_io global)
install(FILES README DESTINATION doc)
