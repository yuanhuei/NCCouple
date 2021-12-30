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
set(VTK_BIN_DIR "./tools/VTK-9.0.1/bin")
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
set(SOURCE_FILE ${SOURCE_FILE1} ${SOURCE_FILE2} ${SOURCE_FILE3} ${SOURCE_FILE4} ${SOURCE_FILE5})
set(INCLUDE_DIR 
	./NCCoouple 
	./tools
	./tools/VTK-9.0.1/include/vtk-9.0
	)

if(WIN32)
	set(STATIC_LIBRARY_DIR ./tools/winlib)
	set(STATIC_LIBRARY_NAME "libmpfr-4" "libgmp-10")
elseif(UNIX)
	set(STATIC_LIBRARY_DIR )
	set(STATIC_LIBRARY_NAME "pthread")
	set(INCLUDE_DIR ${INCLUDE_DIR})
endif()

if (${CMAKE_VERSION} VERSION_LESS 3.8)
	message("Include directories: ${INCLUDE_DIR}")
	include_directories(${INCLUDE_DIR})
	link_directories(${STATIC_LIBRARY_DIR})
	link_libraries(${STATIC_LIBRARY_NAME})
	add_executable (${PROJECT_NAME} ${SOURCE_FILE})
	add_executable (MYTEST "./test/main.cpp")
else()
	add_executable (${PROJECT_NAME} ${SOURCE_FILE})
	target_include_directories(${PROJECT_NAME} PUBLIC ${INCLUDE_DIR})
	target_link_directories(${PROJECT_NAME} PUBLIC ${STATIC_LIBRARY_DIR})
	target_link_libraries(${PROJECT_NAME} PUBLIC ${STATIC_LIBRARY_NAME})

	target_include_directories(MYTEST PUBLIC ${INCLUDE_DIR})
	target_link_directories(MYTEST PUBLIC ${STATIC_LIBRARY_DIR})
	target_link_libraries(MYTEST PUBLIC ${STATIC_LIBRARY_NAME})
endif()