cmake_minimum_required(VERSION 2.8)
project(ARUCO)
set(OpenCV_DIR "/usr/local/opencv/lib/cmake/opencv4")
find_package(OpenCV REQUIRED)
add_executable(main main.cpp)
target_link_libraries(main ${OpenCV_LIBS})
include_directories(~/usr/local/local/opencv/include/opencv4)
