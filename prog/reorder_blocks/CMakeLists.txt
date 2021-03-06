add_executable(reorder_blocks
  reorder_blocks.f90
)
target_link_libraries(reorder_blocks cfdtools_global)
target_link_libraries(reorder_blocks
  progtools
  xyzq_io
)
install(FILES README DESTINATION doc RENAME reorder_blocks)
install(TARGETS reorder_blocks DESTINATION ${CMAKE_INSTALL_BINDIR})
