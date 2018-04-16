add_library(forsythe
  decomp.f
  fmin.f
  fmin77.f
  quanc8.f
  rkf45.f
  seval.f
  solve.f
  spline.f
  svd.f
  urand.f
  zeroin.f
)
target_link_libraries(forsythe global)
install(FILES README DESTINATION doc)
