cfdtools_add_executable(merge_files
  merge_files.f90
)
target_link_libraries(cfdtools_merge_files
  cfdtools_tecplot_io
  cfdtools_progtools
)
