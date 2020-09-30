
if(Boost_FOUND)
  include_directories(${PROJECT_INCLUDE_DIRECTORIES} ${CMAKE_CURRENT_SOURCE_DIR} ${Boost_INCLUDE_DIR})
  link_directories(${PROJECT_LINK_DIRECTORIES} ${Boost_LIBRARY_DIRS})
else()
  include_directories(${PROJECT_INCLUDE_DIRECTORIES} ${CMAKE_CURRENT_SOURCE_DIR})
  link_directories(${PROJECT_LINK_DIRECTORIES})
endif()


add_executable(cbm_qa cbm_qa.cpp)
add_dependencies(cbm_qa AnalysisTreeQA ${PROJECT_DEPENDENCIES})
target_link_libraries(cbm_qa ${ROOT_LIBRARIES} AnalysisTreeQA AnalysisTreeBase AnalysisTreeInfra)

add_executable(pfs_qa pfs_qa.cpp)
add_dependencies(pfs_qa AnalysisTreeQA ${PROJECT_DEPENDENCIES})
target_link_libraries(pfs_qa ${ROOT_LIBRARIES} AnalysisTreeQA AnalysisTreeBase AnalysisTreeInfra)

add_executable(build_effmap build_effmap.cpp)
add_dependencies(build_effmap AnalysisTreeQA ${PROJECT_DEPENDENCIES})
target_link_libraries(build_effmap ${ROOT_LIBRARIES} AnalysisTreeQA AnalysisTreeBase AnalysisTreeInfra)

add_executable(qa3D qa3D.cpp)
add_dependencies(qa3D AnalysisTreeQA ${PROJECT_DEPENDENCIES})
target_link_libraries(qa3D ${ROOT_LIBRARIES} AnalysisTreeQA AnalysisTreeBase AnalysisTreeInfra)

install (TARGETS cbm_qa RUNTIME DESTINATION bin)
install (TARGETS pfs_qa RUNTIME DESTINATION bin)
install (TARGETS build_effmap RUNTIME DESTINATION bin)
install (TARGETS qa3D RUNTIME DESTINATION bin)