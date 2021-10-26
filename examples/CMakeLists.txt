include(AnalysisTree)

if(Boost_FOUND)
    include_directories(${PROJECT_INCLUDE_DIRECTORIES} ${CMAKE_CURRENT_SOURCE_DIR} ${Boost_INCLUDE_DIR} ${AnalysisTree_BINARY_DIR}/include)
    link_directories(${PROJECT_LINK_DIRECTORIES} ${Boost_LIBRARY_DIRS})
else()
    include_directories(${PROJECT_INCLUDE_DIRECTORIES} ${CMAKE_CURRENT_SOURCE_DIR} ${AnalysisTree_BINARY_DIR}/include)
    link_directories(${PROJECT_LINK_DIRECTORIES})
endif()

if(PROJECT_LINK_DIRECTORIES)
    include_directories(${PROJECT_INCLUDE_DIRECTORIES} ${AnalysisTree_BINARY_DIR}/include)
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

add_executable(cbm_qa cbm_qa.cpp)
add_dependencies(cbm_qa AnalysisTreeQA)
target_link_libraries(cbm_qa ${ROOT_LIBRARIES} AnalysisTreeQA AnalysisTreeBase AnalysisTreeInfra)

install (TARGETS cbm_qa RUNTIME DESTINATION bin)

add_executable(covariances_qa covariances_qa.cpp)
add_dependencies(covariances_qa AnalysisTreeQA)
target_link_libraries(covariances_qa ${ROOT_LIBRARIES} AnalysisTreeQA AnalysisTreeBase AnalysisTreeInfra)

install (TARGETS covariances_qa RUNTIME DESTINATION bin)

