cfdtools_add_executable(extract_blayer_data
  extract_blayer_data.f90
)
target_link_libraries(cfdtools_extract_blayer_data
  cfdtools_tecplot_io
  cfdtools_progtools
)
