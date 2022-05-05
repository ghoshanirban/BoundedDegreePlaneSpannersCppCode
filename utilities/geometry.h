//
// Created by matt on 2/5/22.
//

#ifndef LIBSPANNER_GEOMETRY_H
#define LIBSPANNER_GEOMETRY_H

#include "utilities/constants.h"
#include "./globaltypes.h"

namespace BoundedDegreePlaneSpanners {

    template<typename Point_2>
    inline number_t getDistance(const Point_2 &p, const Point_2 &q) {
        return p == q ? 0 : CGAL::sqrt(CGAL::squared_distance(p, q));
    }

    template<class Points>
    inline number_t getEdgeLength(const Edge &e, const Points &P) {
        return getDistance(P[e.first], P[e.second]);
    }

    template<class Point_2>
    inline number_t getAngle(const Point_2 &p, const Point_2 &q, const Point_2 &r) {
        auto pq = q - p,
                rq = q - r;

        number_t result = atan2(pq.y(), pq.x()) - atan2(rq.y(), rq.x());
        result = fmod(result + 2 * PI, 2 * PI);
        return result;
    }

    template<class Points>
    inline number_t getAngle(const index_t p, const index_t q, const index_t r, const Points &P) {
        return getAngle(P[p], P[q], P[r]);
    }

}
#endif //LIBSPANNER_GEOMETRY_H
