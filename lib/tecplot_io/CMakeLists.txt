cfdtools_add_library(tecplot_io
  tec360_dummy.f90
  tecplot_io.f90
  # tecplot_tri_io.f90
)
target_link_libraries(cfdtools_tecplot_io PRIVATE cfdtools_progtools)
