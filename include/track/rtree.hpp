//
// Created by acsr on 9/6/22.
//

#ifndef RACECARPLANNER_RTREE_HPP
#define RACECARPLANNER_RTREE_HPP

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <casadi/casadi.hpp>
#include <vector>

namespace bg = boost::geometry;
namespace bgi=boost::geometry::index;

namespace acsr {
    class RTree {

        typedef bg::model::point<double, 2, bg::cs::cartesian> point_t;
        typedef std::pair<point_t, unsigned> value_t;
    public:
        RTree() = default;
        RTree(const casadi::DM& dm){
            auto n = dm.rows();
            for(auto i=0;i<n;++i){
                tree.insert(std::make_pair(point_t(double(dm(i,0)),double(dm(i,1))),i));
            }
        }

        RTree(const std::vector<std::vector<double>>& dm){
            auto n = dm.size();
            for(auto i=0;i<n;++i){
                tree.insert(std::make_pair(point_t(dm[i][0],dm[i][1]),i));
            }
        }

        std::tuple<std::vector<double>,unsigned> findNearest(const casadi::DM& dm){
            std::vector<value_t> result_n;
            tree.query(bgi::nearest(point_t(double(dm(0)),double(dm(1))), 1), std::back_inserter(result_n));
            auto result = result_n.front();
            return std::make_tuple(std::vector<double>{result.first.get<0>(),result.first.get<1>()},result.second);
        }

        template<typename T>
        std::tuple<std::vector<double>,unsigned> findNearest(const T& dm){
            std::vector<value_t> result_n;
            tree.query(bgi::nearest(point_t(dm[0],dm[1]), 1), std::back_inserter(result_n));
            auto result = result_n.front();
            return std::make_tuple(std::vector<double>{result.first.get<0>(),result.first.get<1>()},result.second);
        }

    private:
        bgi::rtree< value_t, bgi::quadratic<16> > tree;
    };

}

#endif //RACECARPLANNER_RTREE_HPP
