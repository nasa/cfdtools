add_executable(reorder_segments
  reorder_segments.f90
)
target_link_libraries(reorder_segments global)
install(FILES README DESTINATION doc RENAME reorder_segments)
install(TARGETS reorder_segments DESTINATION ${CMAKE_INSTALL_BINDIR})
