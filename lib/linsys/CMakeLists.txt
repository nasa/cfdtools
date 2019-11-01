add_library(linsys
  bloktr.f
  btr4.f
  cholesky.f
  cholesky_factorization.f
  cholesky_solution.f
  chsolve.f
  colslv.f
  decbt.f
  decbtc.f
  decomp.f90
  decslv.f
  dtdlsq.f
  hdecom.f
  hdesol.f
  hdesolw.f
  hsolve.f
  hsulve.f
  lsqr.f
  lusolve.f90
  qrdiag.f
  safeguarded_lsqr.f90
  solbt.f
  solbtc.f
  solve.f90
  symmlq.f
  trdiag.f
  tricps.f
  trid2r.f
  trid3r.f
  trip.f
)
target_link_libraries(linsys PRIVATE global)
target_include_directories(linsys PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
install(FILES README DESTINATION doc RENAME linsys)
