cmake_minimum_required (VERSION 3.14.0)
project (chalk)

include(FetchContent)

FetchContent_Declare(
  cli11
  GIT_REPOSITORY https://github.com/CLIUtils/CLI11
  GIT_TAG        b9be5b9444772324459989177108a6a65b8b2769 # v2.2.0
)
FetchContent_MakeAvailable(cli11)

find_package (Eigen3 3.3 REQUIRED NO_MODULE)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O3 -march=native -Wall -Wextra -Werror -Wno-deprecated-declarations -Wno-unused-function")
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON) # for int128_t

add_subdirectory(util EXCLUDE_FROM_ALL)

include_directories(src)

# the library
file(GLOB files_cpp "src/chalk/*.cpp")
add_library(chalk STATIC ${files_cpp})
target_compile_features(chalk PUBLIC cxx_std_17)
target_link_libraries(chalk util stdc++ gmp mpfr Eigen3::Eigen readline)

# executables
foreach(bin "factor" "primes" "wieferich" "sprp" "calculator" "ziggurat" "exponential_product" "langevin_scheme")
	add_executable(${bin} src/apps/${bin}.cpp)
	target_compile_features(${bin} PUBLIC cxx_std_17)
	target_link_libraries(${bin} chalk CLI11::CLI11)
endforeach(bin)

# unittests
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)

	FetchContent_Declare(
	  Catch2
	  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
	  GIT_TAG        v3.0.1
	)

	FetchContent_MakeAvailable(Catch2)

	add_compile_options(-fsanitize=address -fno-omit-frame-pointer)
	add_link_options(-fsanitize=address)

	add_executable(tests tests/bch.cpp tests/ddouble.cpp tests/integer.cpp tests/langevin_scheme.cpp tests/mpfr_eigen.cpp tests/numerics.cpp tests/rings.cpp)
	target_compile_features(tests PRIVATE cxx_std_17)
	target_link_libraries(tests PRIVATE Catch2::Catch2WithMain chalk)
	target_compile_options(tests PRIVATE -Wall -Wextra -Werror -g -Wno-parentheses -Wno-stringop-overflow)

endif()
