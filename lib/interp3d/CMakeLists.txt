add_library(interp3d
  bilincut.f
  furthest_point.f90
  inside_outside.f90
  intsec4.f
  intsec5.f
  intsec6.f90
  intsec8.f90
  lcsint3d.f
  line_surface.f90
  peval3.f
  pfit3.f
  plbicut.f
  pleval.f
  plfit.f
  plscrv3d.f
  plscurve.f
  plsfit3d.f
  plsintyz.f
  smoothxyz.f
  table3.f
  tri_line.f90
  trilint.f
  wbintr.f
)
target_link_libraries(interp3d global adt)
target_include_directories(interp3d PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
install(FILES README DESTINATION doc)
