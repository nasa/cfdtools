add_executable(extract_columns
  extract_columns.f90
)
target_link_libraries(extract_columns cfdtools_global)
target_link_libraries(extract_columns
  progtools
)
install(FILES README DESTINATION doc RENAME extract_columns)
install(TARGETS extract_columns DESTINATION ${CMAKE_INSTALL_BINDIR})
