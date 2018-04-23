add_executable(gestimate
  gestimate.f90
)
target_link_libraries(gestimate global linsys)
install(FILES README DESTINATION doc RENAME gestimate)
install(TARGETS gestimate DESTINATION ${CMAKE_INSTALL_BINDIR})
