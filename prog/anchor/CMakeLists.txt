add_executable(anchor
  anchor.f90
)
target_link_libraries(anchor global)
target_link_libraries(anchor
  linsys
  numodules
  optinterp
  progtools
)
install(FILES README DESTINATION doc RENAME anchor)
install(TARGETS anchor DESTINATION ${CMAKE_INSTALL_BINDIR})
