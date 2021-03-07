cfdtools_add_executable(combine
  combine.f
)
target_link_libraries(cfdtools_combine
  cfdtools_interplib
  cfdtools_searchlib
  cfdtools_numodules
)
