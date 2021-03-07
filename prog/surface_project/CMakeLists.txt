add_executable(surface_project
  surface_project.f90
)
target_link_libraries(surface_project
  cfdtools_global
  cfdtools_geomlib
  cfdtools_linsys
  cfdtools_progtools
  cfdtools_tecplot_io
  cfdtools_tecplot_adt
)
install(TARGETS surface_project DESTINATION ${CMAKE_INSTALL_BINDIR})
