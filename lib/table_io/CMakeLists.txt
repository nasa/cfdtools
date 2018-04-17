add_library(table_io
  table_io.f90
  test_table_io.f90
)
target_link_libraries(table_io PRIVATE global)
target_include_directories(table_io PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
install(FILES README DESTINATION doc)
