project(circuit-planner)
cmake_minimum_required(VERSION 2.8)
aux_source_directory(. SRC_LIST)
aux_source_directory(src SRC_LIST)
aux_source_directory(src/utils SRC_LIST)
aux_source_directory(extern SRC_LIST)
aux_source_directory(test SRC_LIST)
add_executable(${PROJECT_NAME} ${SRC_LIST})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Werror=return-type")
target_link_libraries(${PROJECT_NAME} stdc++ m)
