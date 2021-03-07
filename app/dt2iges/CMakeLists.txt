cfdtools_add_executable(dt2iges
  dt2iges.f
)
target_link_libraries(cfdtools_dt2iges
  cfdtools_aa_nurbs
  cfdtools_dt_nurbs
  cfdtools_progtools
)
