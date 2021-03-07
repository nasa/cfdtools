cfdtools_add_executable(reflect_blocks_2d
  reflect_blocks_2d.f90
)
target_link_libraries(cfdtools_reflect_blocks_2d
  cfdtools_progtools
  cfdtools_xyzq_io
)
