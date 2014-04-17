include_directories(${CMAKE_SOURCE_DIR}/debyer)

# Run gengetopt
find_program(GENGETOPT "gengetopt")
if (GENGETOPT) 
  execute_process(COMMAND ${GENGETOPT} -i debyer.ggo -uFILE --no-handle-version
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/debyer )
  execute_process(COMMAND ${GENGETOPT} -i extend_cmd.ggo -uFILE --no-handle-version
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/debyer )
else()
  message(FATAL_ERROR "Gengetopt not found!")
endif()

# Add static library and executables

add_library(libdebyer STATIC debyer.c lineio.cc fileio.cc atomtables.c)
target_link_libraries(libdebyer ${LINKLIBS})

add_executable(debyer "cmdline.c" "main.cc")
target_link_libraries(debyer libdebyer)

add_executable(dbr_symmetry dbr_symmetry.cc)
target_link_libraries(dbr_symmetry libdebyer)

add_executable(dbr_extend cells.cc extend_cmd.c dbr_extend.cc)
target_link_libraries(dbr_extend libdebyer)

add_executable(dbr_bonds dbr_bonds.cc cells.cc)
target_link_libraries(dbr_bonds libdebyer)

add_executable(dbr_conv conv.cc)
target_link_libraries(dbr_conv libdebyer)

add_executable(dbr_sic sic.cc)
target_link_libraries(dbr_sic libdebyer)
