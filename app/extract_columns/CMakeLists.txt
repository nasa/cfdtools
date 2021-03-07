cfdtools_add_executable(extract_columns
  extract_columns.f90
)
target_link_libraries(cfdtools_extract_columns
  cfdtools_progtools
)
