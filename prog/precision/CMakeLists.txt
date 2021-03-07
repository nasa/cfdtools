cfdtools_add_executable(precision
  f8pt3.f
  precision.f
)
target_link_libraries(cfdtools_precision
  cfdtools_optlib
)
