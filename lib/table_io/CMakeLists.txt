cfdtools_add_library(table_io
  table_io.f90
  test_table_io.f90
)
target_link_libraries(cfdtools_table_io PRIVATE cfdtools_progtools)
