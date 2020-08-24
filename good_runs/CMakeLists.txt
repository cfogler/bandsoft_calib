set(INCLUDES ${ROOT_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/include )
include_directories(${INCLUDES})


set(CMAKE_VERBOSE_MAKEFILE ON)

add_executable(counts_charge counts_charge.cpp ../calib_helper.cpp )
target_link_libraries(counts_charge bandhit_tree ${ROOT_LIBRARIES} )
install( TARGETS counts_charge
		DESTINATION ${CMAKE_SOURCE_DIR}/bin )
