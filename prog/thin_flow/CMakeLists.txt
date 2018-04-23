add_executable(thin_flow
  thin_flow.f90
)
target_link_libraries(thin_flow global)
target_link_libraries(thin_flow
  xyzq_io
)
install(FILES README DESTINATION doc RENAME thin_flow)
install(TARGETS thin_flow DESTINATION ${CMAKE_INSTALL_BINDIR})
