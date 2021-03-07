add_executable(reshape
  reshape.f
)
target_link_libraries(reshape cfdtools_global)
target_link_libraries(reshape
  cfdtools_progtools
)
install(FILES README DESTINATION doc RENAME reshape)
install(TARGETS reshape DESTINATION ${CMAKE_INSTALL_BINDIR})
