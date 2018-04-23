add_executable(reshape
  reshape.f
)
target_link_libraries(reshape global)
target_link_libraries(reshape
  progtools
)
install(FILES README DESTINATION doc RENAME reshape)
install(TARGETS reshape DESTINATION ${CMAKE_INSTALL_BINDIR})
