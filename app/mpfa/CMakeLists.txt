cfdtools_add_executable(mpfa
  MEDLI_analysis.f90
)
target_link_libraries(cfdtools_mpfa
  cfdtools_interplib
  cfdtools_searchlib
  cfdtools_numodules
)
