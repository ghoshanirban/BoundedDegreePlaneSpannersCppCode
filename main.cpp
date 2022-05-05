/*
 * make main.cpp more user-friendly
 *
 * large pointsets on GitHub?
 * linux cmake
 * github url in the paper
 */

/*
 * Build: cmake CMakeLists.txt
 * Run: ./libspanner
 */

#include <iostream>

#include "algorithms/SpannerAlgorithms.h"
#include "utilities/pointsetgenerators.h"
#include "measurements/SpannerMeasurements.h"

using namespace BoundedDegreePlaneSpanners;

void runAlgorithm(const unsigned i, std::vector<Point> &P, std::vector<Edge> &E) {
    if (i == 1)       BGS05(P, E);
    else if (i == 2)  LW04(P, E);
    else if (i == 3)  BSX09(P, E);
    else if (i == 4)  BGHP10(P, E);
    else if (i == 5)  KPX10(P, E);
    else if (i == 6)  KX12(P, E);
    else if (i == 7)  BCC12<7>(P, E);
    else if (i == 8)  BCC12<6>(P, E);
    else if (i == 9)  BKPX15(P, E);
    else if (i == 10) KPT17(P, E);
    else if (i == 11) BHS18(P, E);
    else std::cout << "Invalid choice!" << std::endl;
}

int main() {
    // choose a distribution and generate a random pointset
    PointGenerator g;

    unsigned n = 1000; // set the number of points to be generated
    std::vector<Point> P;
    g.uni_square(n, 100, P);
    /* Other available generators:
        g.uni_disk(n, 100, P);
        g.normal_clustered(10, n / 10, P);
        g.normal_clustered(1, n, P);
        g.grid_contiguous(n, P);
        g.grid_random(n, P);
        g.annulus(n, P);
        g.galaxy(n, P);
    */

    // Points can also be loaded from files
    //g.loadFromFile("../realworldpointsets/Burma.xy",P);

    std::vector<Edge> E; // every edge is represented using a pair of point-ids which are its endpoints

    CGAL::Real_timer clock;
    clock.start();
    runAlgorithm(2, P, E); // choose an algorithm
    clock.stop();

    std::cout << "Time taken: " << clock.time() << " sec." << std::endl;
    std::cout << "|P| = " << P.size() << std::endl;
    std::cout << "|E| = " << E.size() << std::endl;
    std::cout << "Degree: " << degree(E.begin(), E.end()) << std::endl;
    std::cout << "Avg. degree: " << avgDegree(E.begin(), E.end()) << std::endl;
    std::cout << "Stretch Factor (appx.): " << appxStretchFactor(P.begin(), P.end(),
                                                                 E.begin(), E.end()) << std::endl;
    std::cout << "Stretch Factor (exact.): " << exactStretchFactor(P.begin(), P.end(),
                                                                   E.begin(), E.end()) << std::endl;
    std::cout << "Lightness: " << lightness(P.begin(), P.end(),
                                            E.begin(), E.end()) << std::endl;
    return EXIT_SUCCESS;
}
