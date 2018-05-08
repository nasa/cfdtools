add_library(integrate
  avint.f
  csquad.f
  fsquad.f
  integrate.f
  lcsareas.f
  lcsquad.f
  msquad.f
  quanc8rc.f
  sctrpz.f
  test_2d.f90
  trapezoidal_quadrature_2d.f
)
target_link_libraries(integrate PRIVATE global table_io)
target_include_directories(integrate PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
install(FILES README DESTINATION doc RENAME integrate)
