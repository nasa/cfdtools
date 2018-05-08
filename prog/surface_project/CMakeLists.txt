add_executable(surface_project
  surface_project.f90
)
target_link_libraries(surface_project
  global
  geomlib
  linsys
  progtools
  tecplot_io
  tecplot_adt
)
install(TARGETS surface_project DESTINATION ${CMAKE_INSTALL_BINDIR})
