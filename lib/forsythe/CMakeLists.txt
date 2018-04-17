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
target_include_directories(forsythe PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
install(FILES README DESTINATION doc)
