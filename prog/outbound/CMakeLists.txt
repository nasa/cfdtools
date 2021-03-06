add_executable(outbound
  outbound.f90
  param2d.f
  permute_block.f90
  swap_ijk_conventions.f90
)
target_link_libraries(outbound cfdtools_global)
target_link_libraries(outbound
  xyzq_io
  interp3d
  interplib
  numodules
  progtools
  adt
  interp2d
  linsys
  optlib
  gridlib
  geomlib
  searchlib
)
install(FILES README DESTINATION doc RENAME outbound)
install(TARGETS outbound DESTINATION ${CMAKE_INSTALL_BINDIR})
