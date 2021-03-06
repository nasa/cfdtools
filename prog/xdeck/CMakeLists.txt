add_executable(xdeck
  xdeck.f
)
target_link_libraries(xdeck cfdtools_global)
install(FILES README DESTINATION doc RENAME xdeck)
install(TARGETS xdeck DESTINATION ${CMAKE_INSTALL_BINDIR})
