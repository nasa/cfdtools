cfdtools_add_executable(sorter
  sorter.f
)
target_link_libraries(cfdtools_sorter
  cfdtools_searchlib
  cfdtools_progtools
)
