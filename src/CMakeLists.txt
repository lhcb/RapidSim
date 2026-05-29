file(GLOB RapidSim_sources ${PROJECT_SOURCE_DIR}/src/*.cc)

add_executable(RapidSim.exe
  ${RapidSim_sources}
  ${PROJECT_SOURCE_DIR}/src/RapidSim.C
)

if(EvtGen_FOUND)
  target_compile_definitions(RapidSim.exe PRIVATE RAPID_EVTGEN)
  target_link_libraries(RapidSim.exe PRIVATE ${ROOT_LIBRARIES} ${EVTGEN} ${EVTGENEXT})

  if(CMAKE_CXX_COMPILER_ID MATCHES "Clang" AND CMAKE_SYSTEM_NAME MATCHES "Darwin")
    target_link_options(RapidSim.exe PRIVATE "LINKER:-undefined,dynamic_lookup")
  else()
    target_link_options(RapidSim.exe PRIVATE "LINKER:--unresolved-symbols=ignore-in-shared-libs")
  endif()
else()
  target_link_libraries(RapidSim.exe PRIVATE ${ROOT_LIBRARIES})
endif()

install(TARGETS RapidSim.exe DESTINATION ${CMAKE_INSTALL_BINDIR} RUNTIME DESTINATION bin)