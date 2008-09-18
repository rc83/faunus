include_directories(${CMAKE_SOURCE_DIR}/include)

#----- object files for libfaunus -----
set(objs inputfile io titrate point physconst slump container
    potentials hardsphere group particles analysis species) 
set_source_files_properties(${objs} PROPERTIES LANGUAGE CXX)

#----- locate openbabel ------
if (OPENBABEL_FOUND)
  set(objs ${objs} iobabel)
  include_directories(${OPENBABEL_INCLUDE_DIR})
  set(LINKLIBS ${LINKLIBS} ${OPENBABEL_LIBRARIES})
endif (OPENBABEL_FOUND)

#----- locate gromacs -----
if (FAUNUS_GROMACS AND GMX_FOUND)
  include_directories(${GMX_INCLUDE_DIR})
  if (NOT GMX_MOTIF)
    add_definitions(-DGROMACS)
    set(LINKLIBS ${LINKLIBS} ${GMX_LIBRARIES})
  endif (NOT GMX_MOTIF)
  if (GMX_MOTIF AND X11_FOUND AND MOTIF_FOUND)
    add_definitions(-DGROMACS)
    set(LINKLIBS ${LINKLIBS} ${GMX_LIBRARIES} ${X11_Xt_LIB} ${MOTIF_LIBRARIES})
  endif (GMX_MOTIF AND X11_FOUND AND MOTIF_FOUND)
endif(FAUNUS_GROMACS AND GMX_FOUND)

#----- static faunus library -----
add_library (libfaunus.a STATIC EXCLUDE_FROM_ALL ${objs})
set_target_properties(libfaunus.a PROPERTIES OUTPUT_NAME faunus)
target_link_libraries(libfaunus.a ${LINKLIBS})

#----- dynamic faunus library -----
add_library (libfaunus SHARED ${objs})
set_target_properties(libfaunus PROPERTIES OUTPUT_NAME faunus)
target_link_libraries(libfaunus ${LINKLIBS})
install(TARGETS libfaunus
        LIBRARY DESTINATION lib)

#----- header install target -----
install(DIRECTORY ${CMAKE_SOURCE_DIR}/include/faunus
        DESTINATION include
        FILES_MATCHING PATTERN "*.h"
        PATTERN ".svn" EXCLUDE
        PATTERN "CMakeFiles" EXCLUDE)
