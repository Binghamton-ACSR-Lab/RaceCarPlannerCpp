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

    template<class valid_checker_t>
    class PathPreprocessor {

    public:
        PathPreprocessor(std::shared_ptr<valid_checker_t> valid_checker_ptr, DM &pts, const param_t &vehicle_params,
                         double min_edge_margin = 0.1,double max_edge_margin = 2.0):
                         valid_checker_ptr_(valid_checker_ptr),min_edge_margin_(min_edge_margin),max_edge_margin_(max_edge_margin){
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


        PathPreprocessor(std::shared_ptr<valid_checker_t> valid_checker_ptr, const std::vector<std::vector<double>> &pts, const param_t &vehicle_params,
                         double min_edge_margin = 0.1,double max_edge_margin = 2.0) {
            DM waypoints = DM(pts);
            PathPreprocessor(valid_checker_ptr, waypoints, vehicle_params, min_edge_margin, max_edge_margin);
        }

        std::tuple<bool,std::vector<std::vector<double>>> optimize(std::shared_ptr<Path> path_ptr, int resolution=100,double ds_holder = 0.5){
            auto ds = path_ptr->get_max_length()/resolution;
            if(ds>ds_holder){
                ds=ds_holder;
                resolution = path_ptr->get_max_length()/ds_holder;
            }
            auto s = DM::linspace(0,path_ptr->get_max_length(),resolution+1);
            s = s.T();
            auto t = path_ptr->s_to_t_lookup(s)[0];
            auto zero_n = DM::zeros(1,resolution+1);

            //check centerline
            auto center_line = path_ptr->f_tn_to_xy(DMVector{t,zero_n})[0];
            for(auto i=0;i<resolution+1;++i){
                if(!valid_checker_ptr_->valid(double(center_line(0,i)),double(center_line(1,i))))
                    return std::tuple(false,std::vector<std::vector<double>>{});
            }

            //auto n_array = DM::linspace(-max_edge_margin_,max_edge_margin_,2*steps+1).T();

            //auto inner_n_array2 = DM::linspace(max_edge_margin_+max_edge_margin_/steps,2*max_edge_margin_,steps);
            //auto outer_n_array2 = DM::linspace(-max_edge_margin_-max_edge_margin_/steps,-2*max_edge_margin_,steps);
            //auto inner_n_array1 = DM::linspace(0,max_edge_margin_,steps+1);
            //auto outer_n_array1 = DM::linspace(0,-max_edge_margin_,steps+1);

            auto n_array = DM::linspace(-2*max_edge_margin_,2*max_edge_margin_,4*steps+1).T();
            auto ones = DM::ones(1,4*steps+1).T();

            std::vector<size_t> index_vect(resolution+1);
            std::iota(index_vect.begin(),index_vect.end(),0);


            std::vector<double> inner_boundary_vec(resolution+1,2*max_edge_margin_);
            std::vector<double> outer_boundary_vec(resolution+1,-2*max_edge_margin_);

            std::for_each(std::execution::par,index_vect.begin(),index_vect.end(),[&](size_t idx){
                auto t_array = t(0,idx)*ones;
                //auto inner_pts = path_ptr->f_tn_to_xy(DMVector{t_array,inner_n_array1})[0];
                //auto outer_pts = path_ptr->f_tn_to_xy(DMVector{t_array,outer_n_array1})[0];
                auto pts = path_ptr->f_tn_to_xy(DMVector{t_array,n_array})[0];
                for(auto i=1;i<=2*steps;++i){
                    if(!valid_checker_ptr_->valid(double(pts(0,2*steps+i)),double(pts(1,2*steps+i)))){
                        inner_boundary_vec[idx]=double(n_array(0,2*steps+i-1));
                        break;
                    }
                }

                for(auto i=-1;i>=-2*steps;--i){
                    if(!valid_checker_ptr_->valid(double(pts(0,2*steps+i)),double(pts(1,2*steps+i)))){
                        outer_boundary_vec[idx]=double(n_array(0,2*steps+i+1));
                        break;
                    }
                }
            });

            auto min_inner_dist = *std::min_element(inner_boundary_vec.begin(),inner_boundary_vec.end());
            auto max_outer_dist = *std::max_element(outer_boundary_vec.begin(),outer_boundary_vec.end());
            if(min_inner_dist>=max_edge_margin_ && max_outer_dist<=-max_edge_margin_)
                return std::tuple(true,std::vector<std::vector<double>>{});

            int it = 0;
            std::vector<std::vector<double>> new_pts;
            new_pts.push_back({0,0});

            for(auto i=1;i<resolution;++i){
                if(inner_boundary_vec[i]<min_edge_margin_ && outer_boundary_vec[i]>-min_edge_margin_){
                    return std::tuple(false,std::vector<std::vector<double>>{});
                }

                if(inner_boundary_vec[i]>=max_edge_margin_ && outer_boundary_vec[i]<=-max_edge_margin_){
                    ++it;
                    if(it>=4){
                        new_pts.push_back({double(t(0,i)),0});
                        it=0;
                    }
                }else if(inner_boundary_vec[i]>=max_edge_margin_ && outer_boundary_vec[i]>-max_edge_margin_){
                    if(inner_boundary_vec[i]-outer_boundary_vec[i]>=2*max_edge_margin_){
                        new_pts.push_back({double(t(0,i)),outer_boundary_vec[i]+max_edge_margin_});
                    }else{
                        new_pts.push_back({double(t(0,i)),(outer_boundary_vec[i]+inner_boundary_vec[i])/2});
                    }
                    it=0;
                }else if(inner_boundary_vec[i]<max_edge_margin_ && outer_boundary_vec[i]<=-max_edge_margin_){
                    if(inner_boundary_vec[i]-outer_boundary_vec[i]>=2*max_edge_margin_){
                        new_pts.push_back({double(t(0,i)),inner_boundary_vec[i]-max_edge_margin_});
                    }else{
                        new_pts.push_back({double(t(0,i)),(outer_boundary_vec[i]+inner_boundary_vec[i])/2});
                    }
                    it=0;
                }else{
                    new_pts.push_back({double(t(0,i)),(outer_boundary_vec[i]+inner_boundary_vec[i])/2});
                    it=0;
                }
            }
            new_pts.push_back({double(t(0,-1)),0});
            return std::tuple(true,new_pts);
        }


        DM process(DM &pts) {
            auto total_waypoints = waypoints_.columns();
            for(auto idx=2;idx<total_waypoints-1;++idx) {
                auto idx_vec = split(idx);
                DM new_pts(2,idx+1);

                for(auto i=0;i<idx+1;++i){
                    new_pts(Slice(),i) = waypoints_(Slice(),idx_vec[i]);
                }

                auto path_ptr = std::make_shared<Path>(new_pts,7,50);

                for(auto i=0;i<5;++i){
                    auto [success,new_path_vec] = optimize(path_ptr);
                    if(!success)break;
                    if(success && !new_path_vec.empty()){
                        path_ptr_ = path_ptr;
                        return true;
                    }else{
                        auto new_path = DM(new_path_vec);
                        path_ptr = std::make_shared<Path>(new_path,7,50);
                    }
                }

            }

        }




    private:
        const int steps = 10;
        int segments_;
        double max_edge_margin_;
        double min_edge_margin_;
        DM waypoints_;
        DM dist_vec_;
        std::shared_ptr<valid_checker_t> valid_checker_ptr_;
        std::shared_ptr<Path> path_ptr_;

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
