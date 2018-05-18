add_executable(nozzle_throat_conditions
  nozzle_throat_conditions.f90
)
target_link_libraries(nozzle_throat_conditions
  global
  lewis
  progtools
  integrate1d
  interplib
  searchlib
  numodules
  optlib
  linsys
)
install(FILES README DESTINATION doc RENAME nozzle_throat_conditions)
install(TARGETS nozzle_throat_conditions DESTINATION ${CMAKE_INSTALL_BINDIR})
