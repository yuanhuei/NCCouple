# CMakeList.txt: 顶层 CMake 项目文件，在此处执行全局配置
# 并包含子项目。
#
cmake_minimum_required (VERSION 3.8)

# set the project name and version
set(PROJECT_NAME NCCouple)
project (${PROJECT_NAME} VERSION 1.0)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

aux_source_directory(./NCCouple SOURCE_FILE)
set(INCLUDE_DIR 
	./NCCoouple 
	./tools)

if(WIN32)
	set(STATIC_LIBRARY_DIR ./tools/winlib)
	set(STATIC_LIBRARY_NAME "libmpfr-4" "libgmp-10")
elseif(UNIX)
	set(STATIC_LIBRARY_DIR)
	set(STATIC_LIBRARY_NAME "pthread" "gmp" "mpfr")
endif()

configure_file (
  "./NCCouple/CFDCELLSCoarse.txt"
  "${CMAKE_CURRENT_BINARY_DIR}/CFDCELLSCoarse.txt" 
  COPYONLY
)
configure_file (
  "./NCCouple/pin_c1.apl"
  "${CMAKE_CURRENT_BINARY_DIR}/pin_c1.apl" 
  COPYONLY
)

# 将源代码添加到此项目的可执行文件。
add_executable (${PROJECT_NAME} ${SOURCE_FILE})

target_include_directories(${PROJECT_NAME} PUBLIC ${INCLUDE_DIR})
target_link_directories(${PROJECT_NAME} PUBLIC ${STATIC_LIBRARY_DIR})
target_link_libraries(${PROJECT_NAME} PUBLIC ${STATIC_LIBRARY_NAME})

add_executable (MYTEST "./test/main.cpp")
target_include_directories(MYTEST PUBLIC ${INCLUDE_DIR})
target_link_directories(MYTEST PUBLIC ${STATIC_LIBRARY_DIR})
target_link_libraries(MYTEST PUBLIC ${STATIC_LIBRARY_NAME})