cfdtools_add_executable(extract_functions
  extract_functions.f90
)
target_link_libraries(cfdtools_extract_functions
  cfdtools_progtools
  cfdtools_gridlib
  cfdtools_xyzq_io
)
