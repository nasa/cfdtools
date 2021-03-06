add_executable(combine_blocks
  combine_blocks.f90
)
target_link_libraries(combine_blocks cfdtools_global)
target_link_libraries(combine_blocks
  adt
  linsys
  geomlib
  xyzq_io
)
install(FILES README DESTINATION doc RENAME combine_blocks)
install(TARGETS combine_blocks DESTINATION ${CMAKE_INSTALL_BINDIR})
