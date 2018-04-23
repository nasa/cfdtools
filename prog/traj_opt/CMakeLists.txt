add_executable(traj_opt
  dummies.f
  numerics.f
  second.f
  traj_opt.f
)
target_link_libraries(traj_opt global)
install(FILES README DESTINATION doc RENAME traj_opt)
install(TARGETS traj_opt DESTINATION ${CMAKE_INSTALL_BINDIR})
