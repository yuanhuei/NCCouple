cmake_minimum_required (VERSION 3.0)

# set the project name and version
set(PROJECT_NAME NCCouple)
project (${PROJECT_NAME})

#设置gdb调试
SET(CMAKE_BUILD_TYPE "Debug")  # 定义编译类型
SET(CMAKE_CXX_FLAGS_DEBUG "$ENV{CXXFLAGS} -O0 -Wall -g2 -ggdb") # 定义Debug编译参数
#SET(CMAKE_CXX_FLAGS_RELEASE "$ENV{CXXFLAGS} -O3 -Wall") # 定义Release编译参数

# specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)


aux_source_directory(./NCCouple SOURCE_FILE1)
aux_source_directory(./NCCouple/MHT_common SOURCE_FILE2)
aux_source_directory(./NCCouple/MHT_polyhedron SOURCE_FILE3)
aux_source_directory(./NCCouple/MHT_mesh SOURCE_FILE4)
aux_source_directory(./NCCouple/MHT_field SOURCE_FILE5)
aux_source_directory(./NCCouple/MHT_IO SOURCE_FILE6)
set(SOURCE_FILE ${SOURCE_FILE1} ${SOURCE_FILE2} ${SOURCE_FILE3} ${SOURCE_FILE4} ${SOURCE_FILE5} ${SOURCE_FILE6})
set(INCLUDE_DIR 
	./NCCoouple 
	./tools
	./tools/mpich
	/usr/code/include
	)

if(WIN32)
	#set(Boost_INCLUDE_DIR ${PROJECT_SOURCE_DIR}./tools/boost)
	set(Boost_DIR ${PROJECT_SOURCE_DIR}./tools/boost/lib64-msvc-14.2/lib/cmake/Boost-1.72.0)
	#set(Boost_USE_STATIC_LIBS ON)
	#set(Boost_USE_STATIC_RUNTIME OFF)
	find_package(Boost 1.72.0 COMPONENTS regex)	
	include_directories(include ${Boost_INCLUDE_DIR})
	link_libraries(${Boost_LIBRARIES})
	add_definitions(-D_BasePlatformWinddows_ )
	set(STATIC_LIBRARY_DIR "C:/Program Files (x86)/Microsoft SDKs/MPI/Lib/x64" )
	set(STATIC_LIBRARY_NAME "msmpi")
elseif(UNIX)
	#设置MPI的编译器,不用MPI的话需要去掉
	set (CMAKE_CXX_COMPILER "/opt/hpc/software/mpi/hpcx/v2.4.1/gcc-7.3.1/bin/mpic++")
	add_definitions(-D_BasePlatformLinux_ )
#	set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
	set(STATIC_LIBRARY_DIR "/usr/code/lib")
	set(STATIC_LIBRARY_NAME "pthread" "boost_regex.a") 
	set(INCLUDE_DIR ${INCLUDE_DIR})

endif()

if (${CMAKE_VERSION} VERSION_LESS 3.13)
	message("Include directories: ${INCLUDE_DIR}")
	include_directories(${INCLUDE_DIR})
	link_directories(${STATIC_LIBRARY_DIR})
	link_libraries(${STATIC_LIBRARY_NAME})
	add_executable (${PROJECT_NAME} ${SOURCE_FILE})
else()
	add_executable (${PROJECT_NAME} ${SOURCE_FILE})
	target_include_directories(${PROJECT_NAME} PUBLIC ${INCLUDE_DIR})
	target_link_directories(${PROJECT_NAME} PUBLIC ${STATIC_LIBRARY_DIR})
	target_link_libraries(${PROJECT_NAME} PUBLIC ${STATIC_LIBRARY_NAME})
endif()