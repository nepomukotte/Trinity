cmake_minimum_required(VERSION 3.7)

# set the project name
project(test_utilities)


# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# The BOOST libraries with version > 1.70
find_package(Boost 1.70
             REQUIRED
             COMPONENTS
                system )

# The pthread libreries are need 
find_package(Threads REQUIRED)

include_directories(${PROJECT_SOURCE_DIR}/include /usr/include ${PROJECT_SOURCE_DIR}/../general ${PROJECT_SOURCE_DIR}/../random ${PROJECT_SOURCE_DIR}/../type_traits ${PROJECT_SOURCE_DIR}/../units
 ${PROJECT_SOURCE_DIR}/../string ${PROJECT_SOURCE_DIR}/../chrono ${PROJECT_SOURCE_DIR}/../system ${project_SOURCE_DIR}/../networking ${PROJECT_SOURCE_DIR}/../utility 
 ${PROJECT_SOURCE_DIR}/../serialization ${PROJECT_SOURCE_DIR}/../ ${PROJECT_SOURCE_DIR}/../../ )
# add the executable
add_executable(client_cobo "${PROJECT_SOURCE_DIR}/client_cobo.cpp")
add_executable(client_pdu "${PROJECT_SOURCE_DIR}/client_pdu.cpp")
add_executable(test_random "${PROJECT_SOURCE_DIR}/test_random.cpp")
add_executable(test_typelist "${PROJECT_SOURCE_DIR}/test_typelist.cpp")
add_executable(test_units "${PROJECT_SOURCE_DIR}/test_units.cpp")
add_executable(test_string "${PROJECT_SOURCE_DIR}/test_string.cpp")
add_executable(test_chrono "${PROJECT_SOURCE_DIR}/test_chrono.cpp")
add_executable(test_posix_message_queue "${PROJECT_SOURCE_DIR}/test_posix_message_queue.cpp")
add_executable(test_sync_queue "${PROJECT_SOURCE_DIR}/test_sync_queue.cpp")
add_executable(test_terminal "${PROJECT_SOURCE_DIR}/test_terminal.cpp")
add_executable(test_utility "${PROJECT_SOURCE_DIR}/test_utility.cpp")
add_executable(test_system "${PROJECT_SOURCE_DIR}/test_system.cpp")
add_executable(test_serialization "${PROJECT_SOURCE_DIR}/test_serialization.cpp")
add_executable(test_consumer "${PROJECT_SOURCE_DIR}/test_consumer.cpp")
add_executable(test_consumermanager "${PROJECT_SOURCE_DIR}/test_consumermanager.cpp")
add_executable(test_message "${PROJECT_SOURCE_DIR}/test_message.cpp")
target_link_libraries(test_consumermanager pthread )    
target_link_libraries(test_sync_queue pthread )    
# target_link_libraries(test_utility pthread rt)    
target_link_libraries(test_system PRIVATE Threads::Threads)
if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
	add_executable(master_control "${PROJECT_SOURCE_DIR}/master_control.cpp" "${PROJECT_SOURCE_DIR}/src/state_message_cmds.cpp")
	add_executable(client_siab "${PROJECT_SOURCE_DIR}/client_siab.cpp" "${PROJECT_SOURCE_DIR}/src/commands.cpp")
	add_executable(client_lvps "${PROJECT_SOURCE_DIR}/client_lvps.cpp")
	add_executable(control_software "${PROJECT_SOURCE_DIR}/control_software.cpp")
	add_executable(test_networking "${PROJECT_SOURCE_DIR}/test_networking.cpp")
	add_executable(test_networking2 "${PROJECT_SOURCE_DIR}/test_networking2.cpp")
	add_executable(test_serialport "${PROJECT_SOURCE_DIR}/test_serialport.cpp")
	add_executable(test_serialport2 "${PROJECT_SOURCE_DIR}/test_serialport2.cpp")
	add_executable(test_checker "${PROJECT_SOURCE_DIR}/test_checker.cpp")
	add_executable(test_gcc_comm "${PROJECT_SOURCE_DIR}/test_gcc_comm.cpp")
	add_executable(test_gcc_comm2 "${PROJECT_SOURCE_DIR}/test_gcc_comm2.cpp")
	add_executable(test_tcp_udp_msg "${PROJECT_SOURCE_DIR}/test_tcp_udp_msg.cpp")
	add_executable(test_sequencemanager "${PROJECT_SOURCE_DIR}/test_sequencemanager.cpp")
	add_executable(test_sequencermanager "${PROJECT_SOURCE_DIR}/test_sequencermanager.cpp")
	add_executable(test_process "${PROJECT_SOURCE_DIR}/test_process.cpp")
	add_executable(test_iostream "${PROJECT_SOURCE_DIR}/test_iostream.cpp")
	add_executable(test_monitor "${PROJECT_SOURCE_DIR}/test_monitor.cpp")
# add libraries
	target_link_libraries(master_control boost_system boost_thread boost_program_options pthread )
	target_link_libraries(client_siab boost_system boost_thread boost_program_options pthread )
	target_link_libraries(client_lvps boost_system boost_thread boost_program_options pthread )
	target_link_libraries(control_software boost_system boost_thread boost_program_options pthread rt)
	target_link_libraries(test_networking boost_system boost_thread boost_program_options pthread )    
	target_link_libraries(test_networking2 boost_system boost_thread boost_program_options pthread )    
	target_link_libraries(test_serialport boost_system boost_thread boost_program_options pthread )
	target_link_libraries(test_serialport2 boost_system boost_thread boost_program_options pthread )    
	target_link_libraries(test_gcc_comm boost_system boost_thread boost_program_options pthread)    
	target_link_libraries(test_gcc_comm2 boost_system boost_thread boost_program_options pthread rt)
	target_link_libraries(test_tcp_udp_msg boost_system boost_thread boost_program_options pthread )    
	target_link_libraries(test_sequencemanager boost_system boost_thread boost_program_options pthread )    
	target_link_libraries(test_sequencermanager boost_system boost_thread boost_program_options pthread )    
	target_link_libraries(test_process boost_system boost_thread boost_program_options pthread )    
	target_link_libraries(test_iostream boost_system boost_thread boost_program_options pthread )    
	target_link_libraries(test_monitor boost_system boost_thread boost_program_options pthread )    
endif()

# add libraries
target_link_libraries(client_siab rt)
target_link_libraries(client_cobo rt)
target_link_libraries(client_lvps rt)
target_link_libraries(client_pdu rt)
target_link_libraries(master_control rt)
target_link_libraries(test_posix_message_queue rt)
