cmake_minimum_required(VERSION 3.16)
project(libspanner)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_CXX_FLAGS "-fopenmp -std=c++17 -fexceptions -frounding-math -Wall -Wextra -O3 -Wno-unused-parameter -Wno-unused-but-set-variable")
set(CMAKE_BUILD_TYPE Release)

include_directories(.)

add_executable(libspanner

        algorithms/BCC12.h
        algorithms/BGHP10.h
        algorithms/BGS05.h
        algorithms/BHS18.h
        algorithms/BKPX15.h
        algorithms/BSX09.h
        algorithms/KPT17.h
        algorithms/KPX10.h
        algorithms/KX12.h
        algorithms/LW04.h
        algorithms/types.h
        algorithms/SpannerAlgorithms.h

        delaunay/DelaunayL2.h
        delaunay/DelaunayLinf.h
        delaunay/DelaunayTD.h

        utilities/pointsetgenerators.h
        utilities/ordering.h
        utilities/measure.h
        utilities/mst.h
        utilities/geometry.h
        utilities/constants.h
        utilities/utilities.h
        utilities/globaltypes.h

        measurements/degree.h
        measurements/stretchfactor.h
        measurements/lightness.h
        measurements/SpannerMeasurements.h

        main.cpp )

target_link_libraries(libspanner gmp)
target_link_libraries(libspanner mpfr)
target_link_libraries(libspanner boost_system)
target_link_libraries(libspanner boost_thread)
