cfdtools_add_executable(prepare_local_analysis
  reader.f
  upcase.f
  prepare_local_analysis.f90
)
target_link_libraries(cfdtools_prepare_local_analysis
  cfdtools_progtools
)
