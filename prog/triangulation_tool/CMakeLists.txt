add_executable(triangulation_tool
  triangulation_tool.f90
)
target_link_libraries(triangulation_tool cfdtools_global)
target_link_libraries(triangulation_tool
  cfdtools_triangulation_io
  cfdtools_interp3d
  cfdtools_progtools
  cfdtools_eigenlib
  cfdtools_geomlib
)
install(FILES README DESTINATION doc RENAME triangulation_tool)
install(TARGETS triangulation_tool DESTINATION ${CMAKE_INSTALL_BINDIR})
