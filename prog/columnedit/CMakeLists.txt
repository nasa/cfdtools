add_executable(columnedit
  columnedit.f90
)
target_link_libraries(columnedit cfdtools_global)
target_link_libraries(columnedit
  cfdtools_progtools
)
install(FILES README DESTINATION doc RENAME columnedit)
install(TARGETS columnedit DESTINATION ${CMAKE_INSTALL_BINDIR})
