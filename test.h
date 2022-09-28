//
// Created by xli185 on 9/25/22.
//

#ifndef RACECARPLANNER_TEST_H
#define RACECARPLANNER_TEST_H

#include <boost/geometry.hpp>
#include <iostream>
#include "utility"

void bgi_distance_test(){
    namespace bg = boost::geometry;
    typedef bg::model::point<double, 2, bg::cs::cartesian> point_t;
    typedef bg::model::polygon<point_t> polygon_t;
    polygon_t polygon{{{0.0, 0.0}, {0.0, 5.0}, {5.0, 5.0}, {5.0, 0.0}, {0.0, 0.0}}};

    point_t pt{2,8};


    auto segments = boost::make_iterator_range(bg::segments_begin(polygon), bg::segments_end(polygon));
    for(auto segment : segments){
        //std::cout<<"("<<segment.first->get<0>()<<","<<segment.first->get<1>()<<")\t"
        //                                                                  "("<<segment.second->get<0>()<<","<<segment.second->get<1>()<<")\n";

        auto pt_with_distance = point_to_segment(pt,*segment.first,*segment.second);
        std::cout<<pt_with_distance.distance<<": "<<pt_with_distance.projected_point.get<0>()<<","<<pt_with_distance.projected_point.get<1>()<<'\n';
    }

    auto com_distance = bg::comparable_distance(polygon,pt);
    std::cout<<double(com_distance);


}

void make_planner_test(){

}



#endif //RACECARPLANNER_TEST_H
