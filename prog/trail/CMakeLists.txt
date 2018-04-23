add_executable(trail
  trail.f
)
target_link_libraries(trail global)
install(FILES README DESTINATION doc RENAME trail)
install(TARGETS trail DESTINATION ${CMAKE_INSTALL_BINDIR})
