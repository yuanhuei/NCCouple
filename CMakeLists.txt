# CMakeList.txt: 顶层 CMake 项目文件，在此处执行全局配置
# 并包含子项目。
#
cmake_minimum_required (VERSION 3.8)

# set the project name and version
set(PROJECT_NAME NCCouple)
project (${PROJECT_NAME} VERSION 1.0)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

#aux_source_directory(./NCCouple SOURCE_FILE)
aux_source_directory("./test" SOURCE_FILE)
set(INCLUDE_DIR 
	./NCCoouple 
	./tools)

if(WIN32)
	set(STATIC_LIBRARY_DIR ./tools/winlib)
	set(STATIC_LIBRARY_NAME "libmpfr-4" "libgmp-10")
elseif(UNIX)
	set(STATIC_LIBRARY_DIR /usr/lib/x86_64-linux-gnu)
	set(STATIC_LIBRARY_NAME "gmp gmpxx mpfr")
endif()

# 将源代码添加到此项目的可执行文件。
add_executable (${PROJECT_NAME} ${SOURCE_FILE})

target_include_directories(${PROJECT_NAME} PUBLIC ${INCLUDE_DIR})
if(WIN32)
	target_link_directories(${PROJECT_NAME} PUBLIC ${STATIC_LIBRARY_DIR})
	target_link_libraries(${PROJECT_NAME} PUBLIC ${STATIC_LIBRARY_NAME})
endif()