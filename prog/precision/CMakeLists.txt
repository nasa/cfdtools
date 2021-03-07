add_executable(precision
  f8pt3.f
  precision.f
)
target_link_libraries(precision cfdtools_global)
target_link_libraries(precision
  cfdtools_optlib
)
install(FILES README DESTINATION doc RENAME precision)
install(TARGETS precision DESTINATION ${CMAKE_INSTALL_BINDIR})
