add_executable(reflect_blocks
  reflect_blocks.f90
)
target_link_libraries(reflect_blocks global)
target_link_libraries(reflect_blocks
  progtools
  xyzq_io
)
install(FILES README DESTINATION doc RENAME reflect_blocks)
install(TARGETS reflect_blocks DESTINATION ${CMAKE_INSTALL_BINDIR})
