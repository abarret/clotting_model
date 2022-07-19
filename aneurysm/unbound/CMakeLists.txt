SET(SRCS main.cpp
	 WallSitesMeshMapping.cpp)

ADD_EXECUTABLE(aneurysm-unbound EXCLUDE_FROM_ALL ${SRCS})
SET_TARGET_PROPERTIES(aneurysm-unbound PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/unbound/aneurysm" OUTPUT_NAME main2d)
TARGET_INCLUDE_DIRECTORIES(aneurysm-unbound PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")

TARGET_LINK_LIBRARIES(aneurysm-unbound PRIVATE UnboundClot)
ADD_DEPENDENCIES(aneurysm-unbound UnboundClot)
ADD_DEPENDENCIES(aneurysm aneurysm-unbound)
