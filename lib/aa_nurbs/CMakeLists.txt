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
target_link_libraries(aa_nurbs global)
install(FILES README DESTINATION doc)
