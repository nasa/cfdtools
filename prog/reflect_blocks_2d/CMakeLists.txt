add_executable(reflect_blocks_2d
  reflect_blocks_2d.f90
)
target_link_libraries(reflect_blocks_2d cfdtools_global)
target_link_libraries(reflect_blocks_2d
  cfdtools_progtools
  cfdtools_xyzq_io
)
install(FILES README DESTINATION doc RENAME reflect_blocks_2d)
install(TARGETS reflect_blocks_2d DESTINATION ${CMAKE_INSTALL_BINDIR})
