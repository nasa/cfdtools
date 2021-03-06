add_executable(wingsections
  wingsections.f
)
target_link_libraries(wingsections cfdtools_global)
target_link_libraries(wingsections
  interplib
  integratelib
  progtools
  dt_nurbs
  optlib
  gridlib
  aa_nurbs
)
install(FILES README DESTINATION doc RENAME wingsections)
install(TARGETS wingsections DESTINATION ${CMAKE_INSTALL_BINDIR})
