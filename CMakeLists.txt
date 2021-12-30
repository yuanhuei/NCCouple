cmake_minimum_required (VERSION 3.0)

# set the project name and version
set(PROJECT_NAME NCCouple)
project (${PROJECT_NAME})

# specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

if(WIN32)
set(VTK_DIR "./tools/VTK-9.0.1/lib/cmake/vtk-9.0")
elseif(UNIX)
set(VTK_DIR "./tools/VTK-9.0.1/lib64/cmake/vtk-9.0")
endif()
find_package( VTK REQUIRED ) 
if (NOT VTK_FOUND)
  message("Skipping marchingCube: ${VTK_NOT_FOUND_MESSAGE}")
  return ()
endif()
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
	./tools/VTK-9.0.1/include/vtk-9.0
	)

if(WIN32)
	add_definitions(-D_BasePlatformWinddows_ )
	set(STATIC_LIBRARY_DIR )
	set(STATIC_LIBRARY_NAME ${VTK_LIBRARIES})
elseif(UNIX)
	add_definitions(-D_BasePlatformLinux_ )
	set(STATIC_LIBRARY_DIR )
	set(STATIC_LIBRARY_NAME "pthread"  ${VTK_LIBRARIES})
	set(INCLUDE_DIR ${INCLUDE_DIR})
endif()

if (${CMAKE_VERSION} VERSION_LESS 3.8)
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