cmake_minimum_required (VERSION 3.9.0)
project (nspt)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native -Wall -Wextra -Wshadow -Wnon-virtual-dtor -Werror")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_subdirectory(util)

include_directories(src)
include_directories(util/src)

file(GLOB files_cpp "src/chalk/*.cpp")
add_library(chalk STATIC ${files_cpp})

foreach(bin "test_primes")
	add_executable(${bin} src/tests/${bin}.cpp)
	target_link_libraries(${bin} chalk fmt stdc++)
endforeach(bin)
