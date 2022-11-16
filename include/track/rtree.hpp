//
// Created by acsr on 9/6/22.
//

#ifndef RACECARPLANNER_RTREE_HPP
#define RACECARPLANNER_RTREE_HPP

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <vector>

namespace bg = boost::geometry;
namespace bgi=boost::geometry::index;

namespace acsr {
    class RTree {

        typedef bg::model::point<double, 2, bg::cs::cartesian> point_t;
        typedef std::pair<point_t, double> value_t;
    public:
        RTree() = default;
        RTree(const std::vector<double>& x,const std::vector<double>& y,const std::vector<double>& t){
            auto n = x.size();
            for (auto i = 0; i < n; ++i)
                tree.insert(std::make_pair(point_t(x[i], y[i]), t[i]));
        }


        template<class T>
        std::vector<T> findNearest(const std::vector<T>& vec_x,const std::vector<T>& vec_y){
            std::vector<value_t> result_n;
            std::vector<T> t;
            auto size = vec_x.size();
            for(auto i=0;i<size;++i) {
                tree.query(bgi::nearest(point_t(vec_x[i], vec_y[i]), 1), std::back_inserter(result_n));
            }
            std::transform(result_n.begin(),result_n.end(),t.begin(),[](auto& value){
                return T{value.second};
            });
            return t;
        }

    private:
        bgi::rtree< value_t, bgi::quadratic<16> > tree;
    };

}

#endif //RACECARPLANNER_RTREE_HPP
