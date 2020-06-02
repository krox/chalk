cmake_minimum_required (VERSION 3.10.0)
project (chalk)

include(GoogleTest)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native -Wall -Wextra -Wno-shadow -Wnon-virtual-dtor -Werror -g")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native -Wall -Wextra -Wno-shadow -Wnon-virtual-dtor -Werror -fsanitize=address -fno-omit-frame-pointer -static-libasan -g -fPIE -pie")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON) # for int128_t

add_subdirectory(util)
add_subdirectory(fmt)
add_subdirectory(abseil-cpp)
add_subdirectory(googletest)

include_directories(src)
include_directories(util/src)
include_directories(fmt/include)

# the library
file(GLOB files_cpp "src/chalk/*.cpp")
add_library(chalk STATIC ${files_cpp})
target_link_libraries(chalk fmt stdc++ absl::flat_hash_map absl::inlined_vector gmp mpfr)

# executables
foreach(bin "langevin_scheme" "factor")
	add_executable(${bin} src/apps/${bin}.cpp)
	target_link_libraries(${bin} chalk fmt stdc++ absl::flat_hash_map absl::inlined_vector gmp mpfr)
endforeach(bin)

# tests
enable_testing()
add_subdirectory(src/tests)
