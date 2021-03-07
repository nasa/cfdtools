cfdtools_add_executable(nbyn  nbyn.f)
cfdtools_add_executable(nbyn2 nbyn2.f)
cfdtools_add_executable(nbyn3 nbyn3.f)

target_link_libraries(cfdtools_nbyn  cfdtools_forsythe)
target_link_libraries(cfdtools_nbyn2 cfdtools_linsys)
target_link_libraries(cfdtools_nbyn3 cfdtools_linsys)

