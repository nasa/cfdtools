add_executable(wingsections
  wingsections.f
)
target_link_libraries(wingsections cfdtools_global)
target_link_libraries(wingsections
  cfdtools_interplib
  cfdtools_integratelib
  cfdtools_progtools
  cfdtools_dt_nurbs
  cfdtools_optlib
  cfdtools_gridlib
  cfdtools_aa_nurbs
)
install(FILES README DESTINATION doc RENAME wingsections)
install(TARGETS wingsections DESTINATION ${CMAKE_INSTALL_BINDIR})
