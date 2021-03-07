add_executable(insert_z
  insert_z.f90
)
target_link_libraries(insert_z cfdtools_global cfdtools_tecplot_io cfdtools_progtools)
install(TARGETS insert_z DESTINATION ${CMAKE_INSTALL_BINDIR})
