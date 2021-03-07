add_executable(gestimate
  gestimate.f90
)
target_link_libraries(gestimate cfdtools_global cfdtools_linsys)
install(FILES README DESTINATION doc RENAME gestimate)
install(TARGETS gestimate DESTINATION ${CMAKE_INSTALL_BINDIR})
