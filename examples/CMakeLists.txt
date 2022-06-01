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

add_executable(lambda_qa lambda_qa.cpp)
add_dependencies(lambda_qa AnalysisTreeQA)
target_link_libraries(lambda_qa ${ROOT_LIBRARIES} AnalysisTreeQA AnalysisTreeBase AnalysisTreeInfra)

install (TARGETS lambda_qa RUNTIME DESTINATION bin)

add_executable(tof_purities_qa tof_purities_qa.cpp)
add_dependencies(tof_purities_qa AnalysisTreeQA)
target_link_libraries(tof_purities_qa ${ROOT_LIBRARIES} AnalysisTreeQA AnalysisTreeBase AnalysisTreeInfra)

install (TARGETS tof_purities_qa RUNTIME DESTINATION bin)

add_executable(tof_purities_v0_qa tof_purities_v0_qa.cpp)
add_dependencies(tof_purities_v0_qa AnalysisTreeQA)
target_link_libraries(tof_purities_v0_qa ${ROOT_LIBRARIES} AnalysisTreeQA AnalysisTreeBase AnalysisTreeInfra)

install (TARGETS tof_purities_v0_qa RUNTIME DESTINATION bin)

add_executable(m2_pq m2_pq.cpp)
add_dependencies(m2_pq AnalysisTreeQA)
target_link_libraries(m2_pq ${ROOT_LIBRARIES} AnalysisTreeQA AnalysisTreeBase AnalysisTreeInfra)

install (TARGETS m2_pq RUNTIME DESTINATION bin)

add_executable(v1_qa v1_qa.cpp)
add_dependencies(v1_qa AnalysisTreeQA)
target_link_libraries(v1_qa ${ROOT_LIBRARIES} AnalysisTreeQA AnalysisTreeBase AnalysisTreeInfra)

install (TARGETS v1_qa RUNTIME DESTINATION bin)

