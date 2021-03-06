add_executable(makelist
  makelist.f
)
target_link_libraries(makelist cfdtools_global)
install(FILES README DESTINATION doc RENAME makelist)
install(TARGETS makelist DESTINATION ${CMAKE_INSTALL_BINDIR})
