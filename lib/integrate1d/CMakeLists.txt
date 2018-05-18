add_library(integrate1d
  avint.f
  csquad.f
  fsquad.f
  lcsareas.f
  lcsquad.f
  msquad.f
  quanc8rc.f
  sctrpz.f
  test_2d.f90
  trapezoidal_quadrature_2d.f
)
target_link_libraries(integrate1d PRIVATE global table_io)
target_include_directories(integrate1d PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
install(FILES README DESTINATION doc RENAME integrate1d)
