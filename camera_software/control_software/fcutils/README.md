# fcutils

A collection of general utilities useful to build SPB2 control software and more.

- **algorithm**: 
Same as STL Algorithm. At the moment there is a clamp algorithm … if you do not use C++17 std::clamp.
- **checker**:
A collection of packet checkers and return policies to be used by the serial port, Can bus and TCP/UDP  client/server. Both  
- **chrono**:
Utilities related to the C++11 chrono library. A timer class is present at the moment. 
- **concurrency**:
Utilities related to the C++11 concurrency. A thread joining on deletion as been added.
- **memory**:
Utilities related fo the C++17 polymorpich memory resources (PMR) and other memory mangers or pools.  
- **networking**:
A collection of utilities for the BOOST ASIO library. For example asynchronous TCP and UDP servers. 
- **random**:
Utilities related to the C++11 random generations. A flat and Gaussian generator classes. 
- **serialization**:
Classes to serialize and deserialize plain old data (POD) objects or collection of values.  
- **serialport**:
Utilities related to the BOOST ASIO serial port management. A talk_to_serial port class implements asynchronous reading and synchronous writing to a serial port. 
The class needs two specific types:
	- *Checker*. A functor type implementing the data packet check. It must implement: std::size_t operator()(Data *, std::size), that return the number of bytes to read in the next read operation.
	- *Consumer*. A functor type implementing the consumer of the data. It must implement: void operator()(Data *, std::size), implementing the data usage. 
A little bit of features and checks are needed to be added, but, at this stage, it implement a basic protocol independent class representing a serial communication. A protocol using a data packet terminated with an End Of Line has been added in the test directory. 
- **system**: 
Operating system related facilities. The posix_queue class is saved here, to implement interprocess communication using posix queue.
- **string**:
A collection of utilities needed. A string with added timestamp.
- **terminal**:
A collection of utilities to use terminal colors, management and cursor position using ANSI escape sequences.  
- **type_traits**:
Type traits related utilities not implemented in STL. For example: type_list, is_bool, is_stream_writable, is_stream_readable.
- **units**:
A class implementing MKS group of units and quantities. 
- **utility**:
Am heterogeneous collection of functions and classes. Among them:
	- *extractor and inserter*. Two classes to extract or insert a given number of bits from or into an integer object, given an offset and a mask.
	- *integer*. Type class to obtain integer types given a bit length.
	- *finally*. A class to implement the finally pattern.
	- *overload*. A class to store overload methods for a std::variant visitor.
	- *print_vals*. To print a generalized number of values of different types, comma separated.
	- *operator<< for std::tuple*. To pretty-print a std::tuple to a std::ostream. 
	- *apply_modifier*. To modify contents of a std::tuple applying a function or a functor to them.  
- **test**:
A collection of test programs and examples.  