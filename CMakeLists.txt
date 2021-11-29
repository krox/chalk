cmake_minimum_required (VERSION 3.10.0)
project (chalk)


find_package (Eigen3 3.3 REQUIRED NO_MODULE)

include(GoogleTest)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O3 -march=native -Wall -Wextra -Werror -Wno-deprecated-declarations -Wno-unused-function")
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON) # for int128_t

add_subdirectory(util EXCLUDE_FROM_ALL)
add_subdirectory(abseil-cpp EXCLUDE_FROM_ALL)
add_subdirectory(googletest EXCLUDE_FROM_ALL)

include_directories(src)

# the library
file(GLOB files_cpp "src/chalk/*.cpp")
add_library(chalk STATIC ${files_cpp})
target_compile_features(chalk PUBLIC cxx_std_17)
target_link_libraries(chalk util fmt::fmt stdc++ absl::flat_hash_map absl::inlined_vector gmp mpfr Eigen3::Eigen readline)

# executables
foreach(bin "factor" "primes" "sprp" "calculator" "ziggurat" "exponential_product" "langevin_scheme" "langevin_scheme_1d")
	add_executable(${bin} src/apps/${bin}.cpp)
	target_compile_features(${bin} PUBLIC cxx_std_17)
	target_link_libraries(${bin} chalk)
endforeach(bin)

# tests
enable_testing()
add_subdirectory(src/tests)
