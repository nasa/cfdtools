add_library(table_io
  table_io.f90
  table_io.f90.before-alphanumeric-dataset-option
  test_table_io.f90
)
target_link_libraries(table_io global)
install(FILES README DESTINATION doc)
