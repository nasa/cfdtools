add_library(tecplot_io
  tec360_dummy.f90
  tecplot_io.f90
  tecplot_io.f90.before-bigger-zones
  tecplot_tri_io.f90
  tecplot_tri_io.f90.spare.in.case.of.file.name.mix-up
)
target_link_libraries(tecplot_io global)
install(FILES README DESTINATION doc)
