add_executable(sort_rows
  sort_rows.f90
)
target_link_libraries(sort_rows cfdtools_global)
target_link_libraries(sort_rows
  table_io
  searchlib
  progtools
)
install(FILES README DESTINATION doc RENAME sort_rows)
install(TARGETS sort_rows DESTINATION ${CMAKE_INSTALL_BINDIR})
