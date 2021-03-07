add_executable(dt2iges
  dt2iges.f
)
target_link_libraries(dt2iges
  cfdtools_global
  cfdtools_aa_nurbs
  cfdtools_dt_nurbs
  cfdtools_progtools
)
install(TARGETS dt2iges DESTINATION ${CMAKE_INSTALL_BINDIR})
