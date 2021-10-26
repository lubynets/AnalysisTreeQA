include(AnalysisTree)

set(SOURCES
    EntryConfig.cpp
    Task.cpp
)

string(REPLACE ".cpp" ".hpp" HEADERS "${SOURCES}")
list(APPEND HEADERS "BasicQA.hpp" "CbmCuts.h")

add_library(AnalysisTreeQA SHARED ${SOURCES} G__AnalysisTreeQA.cxx)

if (Boost_FOUND)
  include_directories(${CMAKE_SOURCE_DIR}/src ${PROJECT_INCLUDE_DIRECTORIES} ${Boost_INCLUDE_DIR} ${AnalysisTree_BINARY_DIR}/include)
  link_directories(${PROJECT_LINK_DIRECTORIES} ${Boost_LIBRARY_DIRS})
else ()
  include_directories(${CMAKE_SOURCE_DIR}/src ${PROJECT_INCLUDE_DIRECTORIES} ${AnalysisTree_BINARY_DIR}/include)
  link_directories(${PROJECT_LINK_DIRECTORIES})
endif ()

#add_dependencies(AnalysisTreeQA ${PROJECT_DEPENDENCIES})

ROOT_GENERATE_DICTIONARY(G__AnalysisTreeQA
        ${HEADERS}
        LINKDEF AnalysisTreeQALinkDef.h
        OPTIONS -I${CMAKE_BINARY_DIR}/include
        OPTIONS -I${AnalysisTree_BINARY_DIR}/include
        )
target_link_libraries(AnalysisTreeQA
        PUBLIC
        AnalysisTreeBase
        AnalysisTreeInfra
        )
target_include_directories(AnalysisTreeQA
        PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        )

install(TARGETS AnalysisTreeQA EXPORT AnalysisTreeQATargets
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        RUNTIME DESTINATION bin
        INCLUDES DESTINATION include
        )

install(
        FILES
        ${HEADERS}
        DESTINATION
        include/AnalysisTreeQA
        COMPONENT
        Devel
)

set(PCM_FILE_NAME libAnalysisTreeQA)

install(
        FILES
        "${CMAKE_CURRENT_BINARY_DIR}/${PCM_FILE_NAME}_rdict.pcm"
        DESTINATION
        lib
        OPTIONAL
)

install(
        FILES
        "${CMAKE_CURRENT_BINARY_DIR}/${PCM_FILE_NAME}.rootmap"
        DESTINATION
        lib
        OPTIONAL
)
