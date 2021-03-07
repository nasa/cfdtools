cfdtools_add_executable(deconstruct
  deconstruct.f90
)
target_link_libraries(cfdtools_deconstruct
  cfdtools_xyzq_io
)
