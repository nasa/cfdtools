cfdtools_add_executable(multicut
  multicut.f
)
target_link_libraries(cfdtools_multicut
  cfdtools_numodules
  cfdtools_ugridlib
)
