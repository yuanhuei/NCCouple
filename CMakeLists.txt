cmake_minimum_required (VERSION 3.0)

# set the project name and version
set(PROJECT_NAME NCCouple)
project (${PROJECT_NAME})

# specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

aux_source_directory(./NCCouple SOURCE_FILE1)
aux_source_directory(./NCCouple/MHT_common SOURCE_FILE2)
set(SOURCE_FILE ${SOURCE_FILE1} ${SOURCE_FILE2})
set(INCLUDE_DIR 
	./NCCoouple 
	./tools)

if(WIN32)
	set(STATIC_LIBRARY_DIR ./tools/winlib)
	set(STATIC_LIBRARY_NAME "libmpfr-4" "libgmp-10")
elseif(UNIX)
	set(STATIC_LIBRARY_DIR /usr/code/lib)
	set(STATIC_LIBRARY_NAME "pthread" "gmp" "mpfr")
	set(INCLUDE_DIR ${INCLUDE_DIR} /usr/code/include)
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

if (${CMAKE_VERSION} VERSION_LESS 3.8)
	message("Include directories: ${INCLUDE_DIR}")
	include_directories(${INCLUDE_DIR})
	link_directories(${STATIC_LIBRARY_DIR})
	link_libraries(${STATIC_LIBRARY_NAME})
	add_executable (${PROJECT_NAME} ${SOURCE_FILE})
	add_executable (MYTEST "./test/main.cpp")
else()
	add_executable (${PROJECT_NAME} ${SOURCE_FILE})
	add_executable (MYTEST "./test/main.cpp")

	target_include_directories(${PROJECT_NAME} PUBLIC ${INCLUDE_DIR})
	target_link_directories(${PROJECT_NAME} PUBLIC ${STATIC_LIBRARY_DIR})
	target_link_libraries(${PROJECT_NAME} PUBLIC ${STATIC_LIBRARY_NAME})

	target_include_directories(MYTEST PUBLIC ${INCLUDE_DIR})
	target_link_directories(MYTEST PUBLIC ${STATIC_LIBRARY_DIR})
	target_link_libraries(MYTEST PUBLIC ${STATIC_LIBRARY_NAME})
endif()