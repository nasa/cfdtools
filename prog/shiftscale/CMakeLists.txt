add_executable(shiftscale
  shiftscale.f
)
target_link_libraries(shiftscale cfdtools_global)
target_link_libraries(shiftscale
  cfdtools_numodules
  cfdtools_progtools
)
install(FILES README DESTINATION doc RENAME shiftscale)
install(TARGETS shiftscale DESTINATION ${CMAKE_INSTALL_BINDIR})
