add_library(integrate
  avint.f
  csquad.f
  fsquad.f
  integrate.f
  integrate.f.previous
  lcsareas.f
  lcsquad.f
  msquad.f
  quanc8rc.f
  sctrpz.f
  test_2d.f90
  trapezoidal_quadrature_2d.f
)
target_link_libraries(integrate global)
install(FILES README DESTINATION doc)
