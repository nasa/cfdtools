add_executable(shiftscale
  shiftscale.f
)
target_link_libraries(shiftscale global)
target_link_libraries(shiftscale
  numodules
  progtools
)
install(FILES README DESTINATION doc RENAME shiftscale)
install(TARGETS shiftscale DESTINATION ${CMAKE_INSTALL_BINDIR})
