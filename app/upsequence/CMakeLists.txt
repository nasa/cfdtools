cfdtools_add_executable(upsequence
  upsequence.f90
)
target_link_libraries(cfdtools_upsequence
  cfdtools_xyzq_io
)
