set(INCLUDES ${ROOT_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/include )
include_directories(${INCLUDES})


set(CMAKE_VERBOSE_MAKEFILE ON)
add_executable(initial_bar initial_bar.cpp )
target_link_libraries(initial_bar bandhit_tree ${ROOT_LIBRARIES} )
install( TARGETS initial_bar
		DESTINATION ${CMAKE_SOURCE_DIR}/bin )

add_executable(initial_run initial_run.cpp ../calib_helper.cpp )
target_link_libraries(initial_run bandhit_tree ${ROOT_LIBRARIES} )
install( TARGETS initial_run
		DESTINATION ${CMAKE_SOURCE_DIR}/bin )
