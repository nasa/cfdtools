add_executable(extract_blocks
  extract_blocks.f90
)
target_link_libraries(extract_blocks cfdtools_global)
target_link_libraries(extract_blocks
  cfdtools_xyzq_io
)
install(FILES README DESTINATION doc RENAME extract_blocks)
install(TARGETS extract_blocks DESTINATION ${CMAKE_INSTALL_BINDIR})
