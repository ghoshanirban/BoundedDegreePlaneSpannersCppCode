#ifndef LIBSPANNER_DEGREE_H
#define LIBSPANNER_DEGREE_H

#include <algorithm> // swap, max_element, accumulate
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "utilities/utilities.h"
#include "utilities/globaltypes.h"

namespace BoundedDegreePlaneSpanners {

    template<typename RandomAccessIterator>
    size_t degree(RandomAccessIterator edgesBegin, RandomAccessIterator edgesEnd) {
        typedef typename RandomAccessIterator::value_type EdgeType;
        typedef typename EdgeType::first_type VertexType;

        const std::vector<EdgeType> edges(edgesBegin, edgesEnd);

        if(edges.empty()) {
            return 0;
        }

        std::unordered_map<VertexType, std::unordered_set<VertexType>> adj;
        // for each edge
        for (auto e : edges) {
            auto first = adj.begin();
            std::tie(first, std::ignore) = adj.emplace(e.first, std::unordered_set<VertexType>());
            (*first).second.insert(e.second);

            auto second = adj.begin();
            std::tie(second, std::ignore) = adj.emplace(e.second, std::unordered_set<VertexType>());
            (*second).second.insert(e.first);
        }
        auto max_el = std::max_element(adj.begin(), adj.end(), [&](const auto &lhs, const auto &rhs) {
            return lhs.second.size() < rhs.second.size();
        });
//        cout<<"Largest degree vertex="<<max_el->first<<endl;

        return max_el->second.size();
    }

    template<typename RandomAccessIterator>
    BoundedDegreePlaneSpanners::number_t avgDegree(RandomAccessIterator edgesBegin, RandomAccessIterator edgesEnd) {
        typedef typename RandomAccessIterator::value_type EdgeType;
        typedef typename EdgeType::first_type VertexType;

        const std::vector<EdgeType> edges(edgesBegin, edgesEnd);
        if(edges.empty()) {
            return 0;
        }
        std::unordered_map<VertexType, std::unordered_set<VertexType>> adj;
        // for each edge
        for (auto e : edges) {
            auto first = adj.begin();
            std::tie(first, std::ignore) = adj.emplace(e.first, std::unordered_set<VertexType>());
            (*first).second.insert(e.second);

            auto second = adj.begin();
            std::tie(second, std::ignore) = adj.emplace(e.second, std::unordered_set<VertexType>());
            (*second).second.insert(e.first);
        }
        auto avg = std::accumulate(adj.begin(), adj.end(), 0.0, [&](const BoundedDegreePlaneSpanners::number_t &sum, const auto &current) {
            return sum + current.second.size();
        }) / BoundedDegreePlaneSpanners::number_t(adj.size());

        return avg;
    }
}

#endif // LIBSPANNER_DEGREE_H


