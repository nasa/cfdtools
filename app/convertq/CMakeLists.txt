cfdtools_add_executable(convertq
  convertq.f
)
target_link_libraries(cfdtools_convertq
  cfdtools_interplib
  cfdtools_searchlib
)
