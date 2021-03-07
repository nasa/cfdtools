add_executable(nbyn  nbyn.f)
add_executable(nbyn2 nbyn2.f)
add_executable(nbyn3 nbyn3.f)

target_link_libraries(nbyn  cfdtools_global cfdtools_forsythe)
target_link_libraries(nbyn2 cfdtools_global cfdtools_linsys)
target_link_libraries(nbyn3 cfdtools_global cfdtools_linsys)

install(FILES README DESTINATION doc RENAME nbyn)
install(TARGETS nbyn nbyn2 nbyn3 DESTINATION ${CMAKE_INSTALL_BINDIR})
