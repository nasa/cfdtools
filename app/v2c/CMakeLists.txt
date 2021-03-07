cfdtools_add_executable(v2c
  v2c.f90
)
target_link_libraries(cfdtools_v2c
  cfdtools_progtools
  cfdtools_gridlib
  cfdtools_xyzq_io
)
