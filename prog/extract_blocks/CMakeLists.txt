cfdtools_add_executable(extract_blocks
  extract_blocks.f90
)
target_link_libraries(cfdtools_extract_blocks
  cfdtools_xyzq_io
)
