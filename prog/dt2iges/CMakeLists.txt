add_executable(dt2iges
  dt2iges.f
)
target_link_libraries(dt2iges
  global
  aa_nurbs
  dt_nurbs
  progtools
)
install(TARGETS dt2iges DESTINATION ${CMAKE_INSTALL_BINDIR})
