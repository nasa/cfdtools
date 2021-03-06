add_executable(triangulation_tool
  triangulation_tool.f90
)
target_link_libraries(triangulation_tool cfdtools_global)
target_link_libraries(triangulation_tool
  triangulation_io
  interp3d
  progtools
  eigenlib
  geomlib
)
install(FILES README DESTINATION doc RENAME triangulation_tool)
install(TARGETS triangulation_tool DESTINATION ${CMAKE_INSTALL_BINDIR})
