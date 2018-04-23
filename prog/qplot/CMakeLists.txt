add_executable(qplot
  qplot-ugl.f
  qplotsrc.f
)
target_link_libraries(qplot global)
install(FILES README DESTINATION doc RENAME qplot)
install(TARGETS qplot DESTINATION ${CMAKE_INSTALL_BINDIR})
