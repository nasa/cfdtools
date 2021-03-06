add_executable(extract_lines
  extract_lines.f
)
target_link_libraries(extract_lines cfdtools_global)
install(FILES README DESTINATION doc RENAME extract_lines)
install(TARGETS extract_lines DESTINATION ${CMAKE_INSTALL_BINDIR})
