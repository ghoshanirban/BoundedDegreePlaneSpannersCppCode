#ifndef LIBSPANNER_STRETCHFACTOR_H
#define LIBSPANNER_STRETCHFACTOR_H

#include <algorithm> // swap
#include <deque>
#include <map>
#include <optional>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

//#include <omp.h>

#include "../delaunay/DelaunayL2.h" 
#include "utilities/geometry.h"
#include "../utilities/ordering.h"
#include "Dijkstra.h"
#include "utilities/globaltypes.h"
#include "utilities/utilities.h"

namespace BoundedDegreePlaneSpanners {

    template<typename PointIterator, typename EdgeIterator>
    number_t appxStretchFactor(PointIterator pointsBegin,
                                PointIterator pointsEnd,
                                EdgeIterator edgesBegin,
                                EdgeIterator edgesEnd,
                                const size_t numberOfThreads = 4) {
        std::vector<Point> P(pointsBegin, pointsEnd);
        const std::vector<Edge> E(edgesBegin, edgesEnd);

        if (P.empty() || E.empty()) {
            return 0.0;
        }

        const index_t n = P.size();

        std::vector<std::unordered_set<index_t>> G(n),
                DelG(n);


        for (const Edge &e: E) {

            G[e.first].insert(e.second);
            G[e.second].insert(e.first);
        }

        std::vector<Edge> edgesOfDT;
        std::vector<index_t> index;
        spatialSort<K>(P, index);

        //Step 1: Construct Delaunay triangulation
        DelaunayL2 DT; //DelaunayTriangulationSFH DT(P, edgesOfDT);
        /*Add IDs to the vertex handle. IDs are the number associated to the vertex, also maped as an index in handles.
          (i.e. Vertex with the ID of 10 will be in location [10] of handles.)*/
        FaceHandle hint;
        //cout<<"del:";
        for (size_t entry: index) {
            auto vh = DT.insert(P[entry], hint);
            hint = vh->face();
            vh->info() = entry;
        }
        VertexHandle v_inf = DT.infinite_vertex();

        // Convert DT to adjacency list
        for (auto vit = DT.finite_vertices_begin(); vit != DT.finite_vertices_end(); ++vit) {
            VertexCirculator N = DT.incident_vertices(vit),
                    done(N);
            do {
                if (N != v_inf)
                    DelG[vit->info()].insert(N->info());
            } while (++N != done);
        }

        struct bfsElement {
            typedef index_t level_t;
            index_t vertex;
            level_t level;

            bfsElement(const index_t &vertex, const level_t &level)
                    : vertex(vertex), level(level) {}
        };
        typedef std::deque<bfsElement> BfsQueue;

        typedef std::multimap<number_t, index_t> ShortestPathsQueue;
        typedef ShortestPathsQueue::iterator SPQHandle;


        std::vector<number_t> stretchFactorOfG(numberOfThreads, 0.0);
        std::vector<std::pair<index_t, index_t>> worstPairOfG(numberOfThreads);
        std::vector<std::unordered_map<index_t, number_t>> tracker(n);

        //#pragma omp parallel for num_threads(numberOfThreads) default(shared)
        for (index_t u = 0; u < n; u++) {

            // BFS variables //

            size_t lvl = 0;
            BfsQueue bfs;
            bfs.emplace_back(u, lvl);
            std::unordered_set<index_t> frontier;
            std::vector<bool> known(n, false);
            known[u] = true;

            // Dijkstra variables //

            std::unordered_map<index_t, number_t> shortestPathLength(n);
            shortestPathLength[u] = 0.0;

            ShortestPathsQueue open;
            std::unordered_map<index_t, SPQHandle> openHandle(n);
            openHandle[u] = open.emplace(shortestPathLength[u], u);

            // Iterative Deepening Loop
            number_t t_u = 0.0, t_lvl = 0.0;

            do {
                ++lvl;
                t_u = t_lvl;
                t_lvl = 0.0;

                while (!bfs.empty() && bfs.begin()->level < lvl) {
                    auto v = bfs.begin()->vertex;
                    bfs.pop_front();

                    for (auto w: DelG[v]) {
                        if (!known[w]) {
                            bfs.emplace_back(w, lvl);
                            known[w] = true;
                            if (!contains(shortestPathLength, w)) {
                                frontier.insert(w);
                            } else {
                                auto t_w = shortestPathLength[w] / getDistance(P[u], P[w]);
                                t_lvl = CGAL::max(t_lvl, t_w);
                            }
                        }
                    }
                }

                while (!frontier.empty() && !open.empty()) {
                    auto nextShortestPath = open.begin();
                    index_t v = nextShortestPath->second;
                    shortestPathLength[v] = nextShortestPath->first;
                    open.erase(nextShortestPath);
                    openHandle.erase(v);

                    if (contains(frontier, v)) {
                        frontier.erase(v);
                        auto t_v = shortestPathLength[v] / getDistance(P[u], P[v]);
                        t_lvl = CGAL::max(t_lvl, t_v);
                    }

                    for (auto w: G[v]) {
                        if (!contains(shortestPathLength, w)) {
                            number_t newDist = shortestPathLength[v] + getDistance(P[w], P[v]);

                            if (!contains(openHandle, w) || newDist < openHandle[w]->first) {
                                if (contains(openHandle, w)) {
                                    open.erase(openHandle[w]);
                                }
                                openHandle[w] = open.emplace(newDist, w);
                            }
                        }
                    }
                }
            } while (!bfs.empty() && !(t_u > t_lvl));

            // Done searching neighborhood of u

            if (t_u > stretchFactorOfG[omp_get_thread_num()])
                stretchFactorOfG[omp_get_thread_num()] = t_u;
        }

        return *std::max_element(stretchFactorOfG.begin(), stretchFactorOfG.end());
    }

    // EXACT MEASUREMENT
    template<typename VertexIterator, typename EdgeIterator>
    number_t exactStretchFactor(VertexIterator pointsBegin,
                                            VertexIterator pointsEnd,
                                            EdgeIterator edgesBegin,
                                            EdgeIterator edgesEnd,
                                            const size_t numberOfThreads = 4) {
        //typedef typename VertexIterator::value_type Point_2;

        std::vector<Point> P(pointsBegin, pointsEnd);
        const std::vector<Edge> E(edgesBegin, edgesEnd);

        if (P.empty() || E.empty()) {
            return 0.0;
        }

        const index_t n = P.size();
        std::vector<std::unordered_set<index_t> > G(n, std::unordered_set<index_t>()); // adjacency list

        // Create list of vertices, map to their indices, and adjacency list
        for (auto eit = E.begin(); eit != E.end(); ++eit) {
            auto p = eit->first,
                    q = eit->second;

            G[p].insert(q);
            G[q].insert(p);
        }

        number_t t_max = 0.0;

        // calculate euclidean getDistance between all pairs
        // #pragma omp parallel for reduction( max: t_max ) num_threads(numberOfThreads) default( shared ) // <-- UNCOMMENT FOR PARALLEL, ALSO LINE 9
        for (index_t i = 0; i < n; ++i) {
            // Euclidean distances
            std::vector<number_t> D(n, INF);
            for (index_t j = 0; j < n; ++j) {
                D.at(j) =
                        i == j ? 0 : getDistance(P.at(i), P.at(j));
            }
            // Shortest paths
            std::vector<number_t> ShortestPaths(n, INF);
            std::vector<index_t> Parents(n);
            Dijkstra(i, P, G, ShortestPaths, Parents);

            // Divide each shortest path getDistance by the euclidean distance between the vertices.
            for (size_t j = 0; j < n; ++j) {
                ShortestPaths.at(j) = ( // avoid /0
                        i == j ? 0 : ShortestPaths.at(j) / D.at(j)
                );
            }
            // Find max_t
            auto t_local = max_element(
                    begin(ShortestPaths),
                    end(ShortestPaths)
            );
            if (*t_local > t_max) {
                t_max = *t_local;
            }
        }
        // Find the big mac daddy stretchFactor aka big money
        return t_max;
    }

} // spanner

#endif //LIBSPANNER_STRETCHFACTOR_H
