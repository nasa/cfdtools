add_executable(hyper_aero
  hyper_aero.source.f
)
target_link_libraries(hyper_aero cfdtools_global)
install(FILES README DESTINATION doc RENAME hyper_aero)
install(TARGETS hyper_aero DESTINATION ${CMAKE_INSTALL_BINDIR})
