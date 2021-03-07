cfdtools_add_library(optinterp
  optimal_interpolation.f90
)
target_link_libraries(cfdtools_optinterp PRIVATE cfdtools_linsys)
