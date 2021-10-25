SET(tides_CPP
  BppApplication.cpp
  BppApplicationTools.cpp
  BppAttributesTools.cpp
  BppExceptions.cpp
  BppFileTools.cpp
  BppNestedStringTokenizer.cpp
  BppStringTokenizer.cpp
  BppTextTools.cpp
  HashTable.cpp
  Diploid.cpp
  RealDiploid.cpp
  Child.cpp
  GenomicRanges.cpp
  RecombinationMap.cpp
  Simulator.cpp
  Vcf.cpp
  Tools.cpp
  tides.cpp
)

ADD_EXECUTABLE (tides ${tides_CPP})

SET (tides-targets tides)

FOREACH (target ${tides-targets})
  # Link (static or shared)
  IF (BUILD_STATIC)
    SET_TARGET_PROPERTIES(${target} PROPERTIES LINK_SEARCH_END_STATIC TRUE)
    TARGET_LINK_LIBRARIES(${target} ${BPP_LIBS_STATIC})
    TARGET_LINK_LIBRARIES (${target} ${LIBS})
    TARGET_LINK_LIBRARIES (${target} Threads::Threads)
  ELSE (BUILD_STATIC)
    TARGET_LINK_LIBRARIES(${target} ${BPP_LIBS_SHARED})
    TARGET_LINK_LIBRARIES (${target} ${LIBS})
    TARGET_LINK_LIBRARIES (${target} Threads::Threads)
  ENDIF (BUILD_STATIC)
ENDFOREACH (target)

INSTALL(TARGETS ${tides-targets} DESTINATION ${CMAKE_INSTALL_BINDIR})

