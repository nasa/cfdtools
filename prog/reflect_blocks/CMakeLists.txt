cfdtools_add_executable(reflect_blocks
  reflect_blocks.f90
)
target_link_libraries(cfdtools_reflect_blocks
  cfdtools_progtools
  cfdtools_xyzq_io
)
