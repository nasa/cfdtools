add_executable(xline
  xline.f
)
target_link_libraries(xline global)
target_link_libraries(xline
  progtools
)
install(FILES README DESTINATION doc RENAME xline)
install(TARGETS xline DESTINATION ${CMAKE_INSTALL_BINDIR})
