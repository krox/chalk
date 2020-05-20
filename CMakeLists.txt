cmake_minimum_required (VERSION 3.9.0)
project (chalk)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native -Wall -Wextra -Wno-shadow -Wnon-virtual-dtor -Werror -g")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native -Wall -Wextra -Wno-shadow -Wnon-virtual-dtor -Werror -fsanitize=address -fno-omit-frame-pointer -static-libasan -g -fPIE -pie")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON)

add_subdirectory(util)
add_subdirectory(fmt)
add_subdirectory(abseil-cpp)

include_directories(src)
include_directories(util/src)
include_directories(fmt/include)


file(GLOB files_cpp "src/chalk/*.cpp")
add_library(chalk STATIC ${files_cpp})
target_link_libraries(chalk fmt stdc++ absl::flat_hash_map absl::inlined_vector gmp)

foreach(bin "test_primes" "test_rings" "test_covariant" "test_numerics" "test_algebra" "test_integer")
	add_executable(${bin} src/tests/${bin}.cpp)
	target_link_libraries(${bin} chalk fmt stdc++ absl::flat_hash_map absl::inlined_vector gmp)
endforeach(bin)

foreach(bin "langevin_scheme")
	add_executable(${bin} src/apps/${bin}.cpp)
	target_link_libraries(${bin} chalk fmt stdc++ absl::flat_hash_map absl::inlined_vector gmp)
endforeach(bin)
