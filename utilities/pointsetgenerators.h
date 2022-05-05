#ifndef SPANNERS_GENERATORS_H
#define SPANNERS_GENERATORS_H

#include <iostream>
#include <vector>
#include <chrono>
#include <cmath>
#include <random>
#include <string>
#include <unordered_set>
#include <utility>

#include <boost/functional/hash.hpp>

#include <CGAL/point_generators_2.h>
#include <CGAL/random_convex_set_2.h>

#include "utilities.h"

namespace BoundedDegreePlaneSpanners {


    class PointGenerator {

        inline static number_t m_perturbationValue = 0.0001;

    public:
        PointGenerator() : m_randCgal(   std::random_device{}()   ) {}

        template<class Container, class Random>
        void perturb(Container &P, number_t val, Random rand = Random() ) {
            CGAL::perturb_points_2(P.begin(),P.end(),val,val,rand);
        }

        void uni_square(const index_t n, const double sizeOfSquare, std::vector<Point> &P) {

            typedef CGAL::Random_points_in_square_2<Point, CGAL::Creator_uniform_2<number_t, Point> > Point_generator;
            Point_generator g(sizeOfSquare / 2 , m_randCgal);

            std::unordered_set<Point> P_unique;
            while (( n - P_unique.size()) > 0)
                std::copy_n(g, n, inserter(P_unique, P_unique.end()));

            std::copy(P_unique.begin(), P_unique.end(), back_inserter(P));
            perturb(P, m_perturbationValue);
        }

        void galaxy(const index_t n, std::vector<Point> &P, const unsigned numSpokes = 5) {
            // see https://itinerantgames.tumblr.com/post/78592276402/a-2d-procedural-galaxy-with-c
            //srand(seed());
            const double spokeAngle = 2 * PI / numSpokes,
                    armOffsetMax = 0.5,
                    rotationFactor = 5,
                    perturbationValue = 0.02;

            std::unordered_set<Point> P_unique;

            while (P_unique.size() < n) {
                //for(index_t i=0; i<n; ++i) {
                double distance = randFloat();
                distance = pow(distance, 2);

                double angle = randFloat() * 2 * PI;
                double armOffset = randFloat() * armOffsetMax;
                armOffset -= armOffsetMax / 2;
                armOffset *= (1 / distance);

                double squaredArmOffset = pow(armOffset, 2);
                squaredArmOffset *= -1 * int(armOffset < 0);
                armOffset = squaredArmOffset;

                double rotation = distance * rotationFactor;

                angle = ((unsigned) (angle / spokeAngle)) * spokeAngle;
                angle += armOffset + rotation;

                P_unique.emplace(cos(angle) * distance, sin(angle) * distance);
            }

            std::copy(P_unique.begin(), P_unique.end(), back_inserter(P));
            perturb(P, perturbationValue);
        }

        void uni_disk(const index_t n, const double radiusOfDisk, std::vector<Point> &P) {
            typedef CGAL::Random_points_in_disc_2<Point, CGAL::Creator_uniform_2<number_t, Point> > Point_generator;
            Point_generator g(radiusOfDisk , m_randCgal);

            std::unordered_set<Point> P_unique;
            //size_t remaining;
            while (( n - P_unique.size()) > 0)
                std::copy_n(g, n, inserter(P_unique, P_unique.end()));

            std::copy(P_unique.begin(), P_unique.end(), back_inserter(P));
            perturb(P, m_perturbationValue);
        }

        void normal_clustered(const index_t pointsInACuster,
                              const index_t numberOfClusters,
                              std::vector<Point> &P,
                              const number_t xStdDev = 2.0,
                              const number_t yStdDev = 2.0) {
            std::mt19937 rngX(seed());
            std::mt19937 rngY(seed());
            std::default_random_engine generatorX(rngX()), generatorY(rngY());
            std::normal_distribution<double> distributionX(0.0, xStdDev), distributionY(2.0, yStdDev);

            std::mt19937 rngShift(std::random_device{}());

            std::uniform_int_distribution shiftDistribution(0, INT32_MAX);

            index_t shiftX, shiftY;
            std::unordered_set<std::pair<index_t, index_t>, boost::hash<std::pair<index_t, index_t>>> S;

            std::unordered_set<Point> P_unique;

            for (index_t c = 0; c < numberOfClusters; c++) {
                if (c != 0) {
                    shiftX = shiftDistribution(rngShift) % (20 * numberOfClusters);
                    shiftY = shiftDistribution(rngShift) % (20 * numberOfClusters);

                    while (!(S.find(std::make_pair(shiftX, shiftY)) == S.end())) {
                        shiftX = shiftDistribution(rngShift) % (20 * numberOfClusters);
                        shiftY = shiftDistribution(rngShift) % (20 * numberOfClusters);
                    }
                } else
                    shiftX = shiftY = 0;

                S.insert(std::make_pair(shiftX, shiftY));

                for (index_t i = 0; i < pointsInACuster; i++) {
                    double x = distributionX(generatorX) + (double)shiftX;
                    double y = distributionY(generatorY) + (double)shiftY;
                    P_unique.emplace(x, y);
                }
            }

            std::copy(P_unique.begin(), P_unique.end(), back_inserter(P));
            perturb(P, m_perturbationValue);
        }

        void grid_contiguous(const index_t n, std::vector<Point> &P) {
            CGAL::points_on_square_grid_2(ceil(std::sqrt(n)), n, std::back_inserter(P),
                                          CGAL::Creator_uniform_2<number_t, Point>());
            perturb(P, m_perturbationValue);
        }

        void grid_random(const index_t n, std::vector<Point> &P) {
            std::unordered_set<std::pair<int, int>, boost::hash<std::pair<int, int>>> S;
            std::unordered_set<Point> P_unique;

            std::mt19937 rngX(seed());
            std::mt19937 rngY(seed());
            std::uniform_int_distribution xDistribution(0, (int) ceil(0.7 * (double)n)), yDistribution(0, (int) ceil(0.7 * (double)n));

            index_t count = 0;

            while (count < n) {
                int x = xDistribution(rngX), y = yDistribution(rngY);

                if (S.find(std::make_pair(x, y)) == S.end()) {
                    P_unique.emplace(x, y);
                    S.insert(std::make_pair(x, y));
                    count++;
                }
            }

            std::copy(P_unique.begin(), P_unique.end(), std::back_inserter(P));
            perturb(P, m_perturbationValue);
        }

        void annulus(const index_t n, std::vector<Point> &P, const double r2 = 1000, const double r1 = 800) {
            assert(r2 > r1);
            std::unordered_set<Point> P_unique;

            std::default_random_engine generator(seed());
            std::uniform_real_distribution<double> distributionR(r1, r2), distributionT(0, 1);

            for (index_t i = 0; i < n; i++) {
                double t = 2 * M_PI * distributionT(generator);
                double r = distributionR(generator);
                P_unique.emplace(r * cos(t), r * sin(t));
            }

            std::copy(P_unique.begin(), P_unique.end(), back_inserter(P));
            perturb(P, m_perturbationValue);
        }

        void loadFromFile(const std::string &filename, std::vector<Point> &P) {
            std::ifstream in(filename);

            if (!in.is_open()) {
                std::cout << "Error opening file!\n";
                return;
            }

            number_t x, y;
            while (in >> x >> y)
                P.emplace_back(x, y);

            in.close();

            perturb(P, m_perturbationValue);
        }

    private:
        CGAL::Random m_randCgal;

        size_t static seed() {
            return std::random_device{}();
        }

        double static randFloat() {
            return static_cast<double>(std::random_device{}()) / static_cast<double>(RAND_MAX);
        }

        template<class Container>
        void perturb(Container &P, number_t val) {
            perturb(P, val, m_randCgal); // the function from ../utilities.h
        }
    };

} // spanner

#endif //SPANNERS_GENERATORS_H
