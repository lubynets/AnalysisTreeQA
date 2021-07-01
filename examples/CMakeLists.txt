include(AnalysisTree)

if(Boost_FOUND)
    include_directories(${PROJECT_INCLUDE_DIRECTORIES} ${CMAKE_CURRENT_SOURCE_DIR} ${Boost_INCLUDE_DIR})
    link_directories(${PROJECT_LINK_DIRECTORIES} ${Boost_LIBRARY_DIRS})
else()
    include_directories(${PROJECT_INCLUDE_DIRECTORIES} ${CMAKE_CURRENT_SOURCE_DIR})
    link_directories(${PROJECT_LINK_DIRECTORIES})
endif()

if(PROJECT_LINK_DIRECTORIES)
    include_directories(${PROJECT_INCLUDE_DIRECTORIES})
    link_directories(${PROJECT_LINK_DIRECTORIES})
endif()

add_executable(example example.cpp)
add_dependencies(example AnalysisTreeQA)
target_link_libraries(example ${ROOT_LIBRARIES} AnalysisTreeQA AnalysisTreeBase AnalysisTreeInfra)

install (TARGETS example RUNTIME DESTINATION bin)

add_executable(pfs_qa pfs_qa.cpp)
add_dependencies(pfs_qa AnalysisTreeQA)
target_link_libraries(pfs_qa ${ROOT_LIBRARIES} AnalysisTreeQA AnalysisTreeBase AnalysisTreeInfra)

install (TARGETS pfs_qa RUNTIME DESTINATION bin)