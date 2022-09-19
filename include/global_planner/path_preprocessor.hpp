//
// Created by acsr on 9/18/22.
//

#ifndef RACECARPLANNER_PATH_PREPROCESSOR_HPP
#define RACECARPLANNER_PATH_PREPROCESSOR_HPP

#include "casadi/casadi.hpp"
#include <future>
#include "assert.h"
#include "path.hpp"
#include <execution>

namespace acsr {
    using namespace casadi;
    using param_t = std::map<std::string, double>;

    template<class map_t>
    class PathPreprocessor {

    public:
        PathPreprocessor(map_t *map, DM &pts, const param_t &vehicle_params,
                         double min_edge_margin = 0.1,double max_edge_margin = 2.0):
                         map_(map),min_edge_margin_(min_edge_margin),max_edge_margin_(max_edge_margin){
            //format waypoints
            if (pts.rows() != 2) {
                assert(pts.columns() == 2);
                pts = pts.T();
            }
            auto all = Slice();
            waypoints_ = pts;

            auto diff = pts(Slice(), Slice(1, pts.columns())) - pts(Slice(), Slice(0, -1));
            auto dist = DM::sqrt(diff(0, all) * diff(0, all) + diff(1, all) * diff(1, all));
            dist_vec_ = DM::cumsum(dist);

        }


        PathPreprocessor(map_t *map, const std::vector<std::vector<double>> &pts, const param_t &vehicle_params,
                         double min_edge_margin = 0.1,double max_edge_margin = 2.0) {
            DM waypoints = DM(pts);
            PathPreprocessor(map, waypoints, vehicle_params, min_edge_margin, max_edge_margin);
        }

        DM optimize(std::shared_ptr<Path> path_ptr,map_t* map,int resolution=100,double ds_holder = 0.5){
            auto ds = path_ptr->get_max_length()/resolution;
            if(ds>ds_holder){
                ds=ds_holder;
                resolution = path_ptr->get_max_length()/ds_holder;
            }
            auto s = DM::linspace(0,path_ptr->get_max_length(),resolution+1);
            s = s.T();
            auto t = path_ptr->s_to_t_lookup(s)[0];
            std::for_each(std::execution::par,t->begin(),t->end(),[](double value){

            });
        }


        DM process(DM &pts) {
            auto total_waypoints = waypoints_.columns();
            for(auto idx=2;idx<total_waypoints-1;++idx) {
                auto idx_vec = split(idx);
                DM new_pts(2,idx+1);

                for(auto i=0;i<idx+1;++i){
                    new_pts(Slice(),i) = waypoints_(Slice(),idx_vec[i]);
                }

                path_ptr = std::make_shared<Path>(new_pts,7,50);

            }

        }


    private:
        int segments_;
        double max_edge_margin_;
        double min_edge_margin_;
        DM waypoints_;
        DM dist_vec_;
        map_t *map_;
        std::shared_ptr<Path> path_ptr;

        std::vector<int> split(int segment) {
            auto all = Slice();
            if (segment == 1)return std::vector<int>{0,waypoints_.columns()-1};

            std::vector<int> result(segment + 1);
            result[0]=0;
            result[segment] = waypoints_.columns()-1;

            auto segment_dist = double(dist_vec_(0, -1)) / segment;
            int idx = 0;
            for (auto i = 0; i < dist_vec_.columns()-1; ++i) {
                if ((double(dist_vec_(0, i)) <= segment_dist * (idx + 1) &&
                    double(dist_vec_(0, i + 1)) >= segment_dist * (idx + 1)) || dist_vec_.columns()-i-1<=segment-idx) {
                    result[idx+1] = i;
                    ++idx;
                    if (idx == segment - 1)return result;
                }
            }
        }


    };
}


#endif //RACECARPLANNER_PATH_PREPROCESSOR_HPP
