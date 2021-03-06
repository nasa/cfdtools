add_executable(insert_z
  insert_z.f90
)
target_link_libraries(insert_z cfdtools_global tecplot_io progtools)
install(TARGETS insert_z DESTINATION ${CMAKE_INSTALL_BINDIR})
