cfdtools_add_executable(insert_z
  insert_z.f90
)
target_link_libraries(cfdtools_insert_z
  cfdtools_tecplot_io
  cfdtools_progtools
)
