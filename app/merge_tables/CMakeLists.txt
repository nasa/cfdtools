cfdtools_add_executable(merge_tables
  merge_tables.f90
)
target_link_libraries(cfdtools_merge_tables
  cfdtools_progtools
)
