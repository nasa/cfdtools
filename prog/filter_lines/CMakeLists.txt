add_executable(filter_lines
  filter_lines.f90
)
target_link_libraries(filter_lines cfdtools_global)
target_link_libraries(filter_lines
  cfdtools_progtools
)
install(FILES README DESTINATION doc RENAME filter_lines)
install(TARGETS filter_lines DESTINATION ${CMAKE_INSTALL_BINDIR})
