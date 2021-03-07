cfdtools_add_executable(thin_flow
  thin_flow.f90
)
target_link_libraries(cfdtools_thin_flow
  cfdtools_xyzq_io
)
