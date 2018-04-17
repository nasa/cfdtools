add_library(interp2d
  akima.f
  bicubmat.f
  bilinear.f
  bilint.f
  bilintp.f
  bvlsq.f
  facets.f
  hardy.f
  inmesh.f
  intsec2.f
  intsec2d.f90
  intsec2t.f
  intsec3.f
  lcsfit2d.f
  lcsint2d.f
  lint2d.f
  nearpt.f
  pbilint.f
  pbvlsq.f
  plbicube.f
  plincub.f
  plsint2d.f
  plstangent.f
  plxcut.f
  smooth2d.f
  table2.f
  tfif2d.f
  tpspln.f
)
target_link_libraries(interp2d global adt)
target_include_directories(interp2d PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
install(FILES README DESTINATION doc)
