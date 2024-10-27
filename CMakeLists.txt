cmake_minimum_required(VERSION 3.16.0)
project(chalk)

include(cmake/CPM.cmake)

CPMAddPackage(
	NAME util
	GITHUB_REPOSITORY krox/util
	GIT_TAG 2c9b21ea14511475f9c8a4d993c19bebf89971d7
)

CPMAddPackage(
	NAME cli11
	GITHUB_REPOSITORY CLIUtils/CLI11
	GIT_TAG f4d0731cebb123ff0ace712c099dffbcd2c58e5a # v2.4.1
)

CPMAddPackage(
	NAME Catch2
	GITHUB_REPOSITORY catchorg/Catch2
	GIT_TAG 2ab20a0e008845e02bd06248e61ca6e5ad1aba33 # v3.3.1
)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O3 -march=native -Wall -Wextra -Werror")
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON) # for int128_t

include_directories(src)

# the library
file(GLOB files_cpp "src/chalk/*.cpp")
add_library(chalk STATIC ${files_cpp})
target_compile_features(chalk PUBLIC cxx_std_17)
target_link_libraries(chalk util stdc++ gmp mpfr Eigen3::Eigen readline)

# executables
foreach(bin "factor" "primes" "wieferich" "psp" "calculator" "ziggurat" "exponential_product" "langevin_scheme")
	add_executable(${bin} src/apps/${bin}.cpp)
	target_compile_features(${bin} PUBLIC cxx_std_17)
	target_link_libraries(${bin} chalk CLI11::CLI11)
endforeach(bin)

# unittests
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
	add_compile_options(-fsanitize=address -fno-omit-frame-pointer)
	add_link_options(-fsanitize=address)

	add_executable(tests tests/bch.cpp tests/ddouble.cpp tests/integer.cpp tests/langevin_scheme.cpp tests/mpfr_eigen.cpp tests/numerics.cpp tests/rings.cpp)
	target_compile_features(tests PRIVATE cxx_std_17)
	target_link_libraries(tests PRIVATE Catch2::Catch2WithMain chalk)
	target_compile_options(tests PRIVATE -Wall -Wextra -Werror -g)
endif()
