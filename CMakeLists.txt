cmake_minimum_required(VERSION 3.28.0)
project(chalk LANGUAGES CXX C)

option(CHALK_BUILD_APPS "build chalk command-line tools" ${PROJECT_IS_TOP_LEVEL})
option(CHALK_BUILD_TESTS "build chalk unit tests" ${PROJECT_IS_TOP_LEVEL})
set(CHALK_COMPILE_OPTIONS "" CACHE STRING "extra compiler flags applied only to chalk-owned targets")

# Disable scanning for C++20 modules. Decreases build times a lot because currently
# * module dependency scanning introduces a lot of false dependencies
# * does not play well with ccache
set(CMAKE_CXX_SCAN_FOR_MODULES OFF)

# Set C++23 explicitly. Note that 'target_compile_features(... cxx_std_23)' alone is not
# sufficient because it does not propagate upstream to the util library.
set(CMAKE_CXX_STANDARD 23)

# GNU extensions are required for __int128.
set(CMAKE_CXX_EXTENSIONS ON)

include(cmake/CPM.cmake)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)

CPMAddPackage(
	NAME util
	GITHUB_REPOSITORY krox/util
	GIT_TAG 1718e6cc30705ba54c0c8c5f8d0c346565f372fb
	OPTIONS
		"UTIL_BUILD_APPS OFF"
		"UTIL_BUILD_TESTS OFF"
)

if(CHALK_BUILD_APPS)
	CPMAddPackage(
		NAME cli11
		GITHUB_REPOSITORY CLIUtils/CLI11
		GIT_TAG 6c7b07a878ad834957b98d0f9ce1dbe0cb204fc9 # v2.4.2
	)
endif()

set(chalk_lib_sources
	src/chalk/indexed.cpp
	src/chalk/numerics.cpp
	src/chalk/numtheory.cpp
	src/chalk/parser.cpp
	src/chalk/symbolic.cpp
)

add_library(chalk STATIC ${chalk_lib_sources})
target_include_directories(chalk PUBLIC src)
target_compile_features(chalk PUBLIC cxx_std_23)
set_target_properties(chalk PROPERTIES CXX_EXTENSIONS ON)
target_link_libraries(chalk PUBLIC util Eigen3::Eigen gmp mpfr)
target_compile_options(chalk PRIVATE ${CHALK_COMPILE_OPTIONS})

if(CHALK_BUILD_APPS)
	foreach(bin factor primes wieferich psp calculator ziggurat exponential_product langevin_scheme)
		add_executable(${bin} src/apps/${bin}.cpp)
		target_link_libraries(${bin} PRIVATE chalk CLI11::CLI11)
		target_compile_options(${bin} PRIVATE ${CHALK_COMPILE_OPTIONS})
	endforeach()
endif()

if(CHALK_BUILD_TESTS)
	include(CTest)
	if(BUILD_TESTING)
		CPMAddPackage(
			NAME Catch2
			GITHUB_REPOSITORY catchorg/Catch2
			GIT_TAG fa43b77429ba76c462b1898d6cd2f2d7a9416b14 # v3.7.1
			OPTIONS "CATCH_INSTALL_DOCS OFF"
		)
		include(${Catch2_SOURCE_DIR}/extras/Catch.cmake)

		add_executable(tests
			tests/bch.cpp
			tests/ddouble.cpp
			tests/integer.cpp
			tests/langevin_scheme.cpp
			tests/mpfr_eigen.cpp
			tests/numerics.cpp
			tests/rings.cpp
			tests/symbolic.cpp
		)
		target_link_libraries(tests PRIVATE Catch2::Catch2WithMain chalk)
		target_compile_features(tests PRIVATE cxx_std_23)
		target_compile_options(tests PRIVATE ${CHALK_COMPILE_OPTIONS})
		catch_discover_tests(tests)
	endif()
endif()
