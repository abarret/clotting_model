SET(SRCS main.cpp
	 WallSitesMeshMapping.cpp)

ADD_EXECUTABLE(aneurysm EXCLUDE_FROM_ALL ${SRCS})
SET_TARGET_PROPERTIES(aneurysm PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/aneurysm" OUTPUT_NAME main2d)
TARGET_INCLUDE_DIRECTORIES(aneurysm PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")

TARGET_LINK_LIBRARIES(aneurysm PRIVATE Clot)
ADD_DEPENDENCIES(aneurysm Clot)
