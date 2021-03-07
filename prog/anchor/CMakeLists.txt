add_executable(anchor
  anchor.f90
)
target_link_libraries(anchor cfdtools_global)
target_link_libraries(anchor
  cfdtools_linsys
  cfdtools_numodules
  cfdtools_optinterp
  cfdtools_progtools
)
install(FILES README DESTINATION doc RENAME anchor)
install(TARGETS anchor DESTINATION ${CMAKE_INSTALL_BINDIR})
