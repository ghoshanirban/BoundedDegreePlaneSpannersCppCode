## Bounded-degree plane geometric spanners in practice

The construction of bounded-degree plane geometric spanners has been a focus of interest  since 2002 when Bose, Gudmundsson, and Smid proposed the first algorithm to construct such spanners. To date, eleven algorithms have been designed with various trade-offs in degree and stretch factor. We have implemented these sophisticated algorithms in C++ using the CGAL library and experimented with them using large synthetic and real-world pointsets. Our  experiments  have revealed their practical behavior and real-world efficacy. We share the implementations via GitHub for broader uses and future research.

We present a simple practical algorithm, named AppxStretchFactor, that can estimate  stretch factors (obtains  lower bounds on the exact stretch factors) of geometric spanners - a challenging problem for which no practical algorithm is known yet. In our experiments with bounded-degree plane geometric spanners, we find that AppxStretchFactor estimates stretch factors almost precisely. Further, it gives linear runtime performance in practice for the pointset distributions considered in this work, making it much faster than the naive Dijkstra-based algorithm for calculating stretch factors

Pre-print: https://arxiv.org/abs/2205.03204

## Built with

To build the project, use CMake. Before running CMake, make sure CGAL and the following required libraries are installed on your system. Please note that depending on your system, the supplied makefile CMakeLists.txt may need to be edited slightly. We have tested the project on Ubuntu 20.04 LTS. With some slight modifications in the CMakeLists.txt, the project can be built on macOS as well. 

* [GCC](https://gcc.gnu.org/)
* [C++17](https://en.cppreference.com/w/cpp/17)
* [CGAL](https://www.cgal.org/)
* [Boost](https://www.boost.org/)
* [GMP](https://gmplib.org/)
* [MPFR](https://www.mpfr.org/)



## Authors

* [Fred Anderson, University of North Florida, FL, USA](https://github.com/TheDKG)
* [Anirban Ghosh, University of North Florida, FL, USA](https://github.com/ghoshanirban)
* [Matthew Graham, University of North Florida, FL, USA](https://github.com/mgatc)
* [Lucas Mougeot, University of North Florida, FL, USA](https://github.com/lucasfuturist)
* [David Wisnosky, University of North Florida, FL, USA](https://github.com/Wisno33)

## Acknowledgement

Research on this project was supported by the University of North Florida Academic Technology Grant and NSF Award CCF-1947887.
