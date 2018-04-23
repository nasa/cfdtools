add_executable(wingsections
  wingsections.f
)
target_link_libraries(wingsections global)
target_link_libraries(wingsections
  interplib
  integrate
  progtools
  dt_nurbs
  optlib
  gridlib
  aa_nurbs
)
install(FILES README DESTINATION doc RENAME wingsections)
install(TARGETS wingsections DESTINATION ${CMAKE_INSTALL_BINDIR})
