add_executable(combine_blocks_turb
  combine_blocks_turb.f90
)
target_link_libraries(combine_blocks_turb cfdtools_global)
target_link_libraries(combine_blocks_turb
  cfdtools_adt
  cfdtools_linsys
  cfdtools_geomlib
  cfdtools_xyzq_io
)
install(FILES README DESTINATION doc RENAME combine_blocks_turb)
install(TARGETS combine_blocks_turb DESTINATION ${CMAKE_INSTALL_BINDIR})
