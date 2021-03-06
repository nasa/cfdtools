add_library(integratelib
  avint.f
  csquad.f
  fsquad.f
  lcsareas.f
  lcsquad.f
  msquad.f
  quanc8rc.f
  sctrpz.f
  trapezoidal_quadrature_2d.f
  trapezoidal_quadrature_3d.f90
)
target_link_libraries(integratelib PRIVATE cfdtools_global table_io interplib)
target_include_directories(integratelib PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
install(FILES README DESTINATION doc RENAME integratelib)
