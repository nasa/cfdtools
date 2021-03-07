cfdtools_add_executable(sort_rows
  sort_rows.f90
)
target_link_libraries(cfdtools_sort_rows
  cfdtools_table_io
  cfdtools_searchlib
  cfdtools_progtools
)
