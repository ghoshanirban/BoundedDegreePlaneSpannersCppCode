#ifndef LIBSPANNER_LIGHTNESS_H
#define LIBSPANNER_LIGHTNESS_H

#include <tuple>

#include "utilities/geometry.h"
#include "../utilities/mst.h"
#include "utilities/globaltypes.h"
#include "utilities/utilities.h"

namespace BoundedDegreePlaneSpanners {

    template< class VertexIterator, class EdgeIterator>
    number_t weight( VertexIterator pointsBegin,
                     VertexIterator pointsEnd,
                     EdgeIterator edgesBegin,
                     EdgeIterator edgesEnd ) {
        const std::vector<Point> P(pointsBegin,pointsEnd);

        number_t w = 0.0;
        index_t p,q;
        for (auto e = edgesBegin; e != edgesEnd; ++e) {
            std::tie(p,q) = *e;
            w += getDistance(P[p],P[q]);
        }
        return w;
    }

    template< class VertexIterator, class EdgeIterator>
    number_t lightness( VertexIterator pointsBegin,
                           VertexIterator pointsEnd,
                           EdgeIterator edgesBegin,
                           EdgeIterator edgesEnd ) {
        std::vector<Point> P(pointsBegin,pointsEnd);
        std::vector<Edge> E(edgesBegin,edgesEnd);
        std::list<Edge> MST;
        getMST( P.begin(), P.end(), E.begin(), E.end(), back_inserter(MST) );
        number_t weightOfMST = weight(P.begin(), P.end(), MST.begin(), MST.end() ),
                weightOfG   = weight(P.begin(), P.end(), E.begin(), E.end() ),
                lightness = weightOfG / weightOfMST;
        return lightness;
    }

} // namespace spanner

#endif // LIBSPANNER_LIGHTNESS_H


