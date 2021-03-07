cfdtools_add_executable(reorder_blocks
  reorder_blocks.f90
)
target_link_libraries(cfdtools_reorder_blocks
  cfdtools_progtools
  cfdtools_xyzq_io
)
