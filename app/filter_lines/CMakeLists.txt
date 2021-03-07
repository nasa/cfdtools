cfdtools_add_executable(filter_lines
  filter_lines.f90
)
target_link_libraries(cfdtools_filter_lines
  cfdtools_progtools
)
