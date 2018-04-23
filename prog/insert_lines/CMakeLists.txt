add_executable(insert_lines
  insert_lines.f90
)
target_link_libraries(insert_lines global)
install(FILES README DESTINATION doc RENAME insert_lines)
install(TARGETS insert_lines DESTINATION ${CMAKE_INSTALL_BINDIR})
