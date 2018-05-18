add_library(aa_nurbs
  bsaf2xy.f
  bsafgeom.f
  bsafnorm.f
  bsafy4x.f
  bsarcmap.f
  bsarea.f
  bscrvmod.f
  bscrvrd.f
  bseval.f
  bsfit.f
  bsinsert.f
  bszeromn.f
  conarea.f
  dtlsa.f
  dtradc.f
  dtrdwr.f
  dtsepp.f
  lsknots.f
)
target_link_libraries(aa_nurbs PRIVATE
  global
  dt_nurbs
  geomlib
  gridlib
  interplib
  integrate1d
  optlib
)
target_include_directories(aa_nurbs PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
install(FILES README DESTINATION doc RENAME aa_nurbs)
