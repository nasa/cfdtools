add_library(optlib
  cendif.f
  cubic.f
  fmin77.f
  fminrc.f
  fminrc2.f
  hooke_jeeves.f90
  hooke_jeeves.f90.unconstrained
  mna.f
  objeps.f
  qnm_module.f90
  qnm_solution.f90
  qnm_solution3.f90
  qnmdif.f
  qnmdif2.f
  zerorc.f
  zeros.f
)
target_link_libraries(optlib global)
install(FILES README DESTINATION doc)
