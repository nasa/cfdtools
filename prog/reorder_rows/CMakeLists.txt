add_executable(reorder_rows
  reorder_rows.f90
)
target_link_libraries(reorder_rows global)
install(FILES README DESTINATION doc RENAME reorder_rows)
install(TARGETS reorder_rows DESTINATION ${CMAKE_INSTALL_BINDIR})
