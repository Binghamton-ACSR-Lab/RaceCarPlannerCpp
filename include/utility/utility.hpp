//
// Created by lab on 8/31/22.
//

#ifndef RACECARPLANNER_UTILITY_HPP
#define RACECARPLANNER_UTILITY_HPP
#include <boost/geometry.hpp>
#include <casadi/casadi.hpp>
namespace bg=boost::geometry;

template<typename P>
struct DistancePoint {
    double distance;
    P projected_point;
};

template <typename P,typename Strategy = bg::strategy::distance::pythagoras<> >
DistancePoint<P> point_to_segment(P const& p, P const& p1, P const& p2) {
    P v = p2, w = p;
    bg::subtract_point(v, p1);
    bg::subtract_point(w, p1);

    auto const c1 = bg::dot_product(w, v);
    if (c1 <= 0)  return { Strategy::apply(p, p1), p1 };

    auto const c2 = bg::dot_product(v, v);
    if (c2 <= c1) return { Strategy::apply(p, p2), p2 };

    P prj = p1;
    bg::multiply_value(v, c1/c2);
    bg::add_point(prj, v);

    return { Strategy::apply(p, prj), prj };
}

template <typename P,typename Strategy = bg::strategy::distance::pythagoras<> >
DistancePoint<P> point_to_polygon(P const& pt, bg::model::polygon<P> const& polygon) {
    DistancePoint<P> distance_point;
    distance_point.distance = std::numeric_limits<double>::max();
    auto segments = boost::make_iterator_range(bg::segments_begin(polygon), bg::segments_end(polygon));
    for(auto segment : segments){
        auto pt_with_distance = point_to_segment(pt,*segment.first,*segment.second);
        if(pt_with_distance.distance<distance_point.distance)
            distance_point =pt_with_distance;
        //std::cout<<pt_with_distance.distance<<": "<<pt_with_distance.projected_point.get<0>()<<","<<pt_with_distance.projected_point.get<1>()<<'\n';
    }
    return distance_point;
}

namespace acsr{
    template<typename T>
    void printArray(T& arr){
        for(auto& t:arr){
            std::cout<<t<<'\t';
        }
        std::cout<<'\n';
    }
}

#endif //RACECARPLANNER_UTILITY_HPP
