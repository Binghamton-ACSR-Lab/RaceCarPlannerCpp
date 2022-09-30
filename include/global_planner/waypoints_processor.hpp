//
// Created by xli185 on 9/25/22.
//

#ifndef RACECARPLANNER_WAYPOINTS_PROCESSOR_HPP
#define RACECARPLANNER_WAYPOINTS_PROCESSOR_HPP


#include "casadi/casadi.hpp"
#include <future>
#include "assert.h"
#include "path.hpp"
#include <execution>

namespace acsr {
    using namespace casadi;
    using param_t = std::map<std::string, double>;
    typedef bg::model::point<double, 2, bg::cs::cartesian> point_t;

    template<class valid_checker_t>
    class WaypointsProcessor {

    public:
        WaypointsProcessor(std::shared_ptr<valid_checker_t> valid_checker_ptr,
                           DM pts,
                           double min_edge_margin = 1,double max_edge_margin = 8.0)
                           :valid_checker_ptr_(valid_checker_ptr),min_edge_margin_(min_edge_margin),max_edge_margin_(max_edge_margin){
            //format waypoints
            if (pts.rows() != 2) {
                assert(pts.columns() == 2);
                pts = pts.T();
            }
            waypoints_ = pts;
            std::cout<<"Construct path from: \n"<<waypoints_<<std::endl;
            original_path_ptr_=std::make_shared<Path>(pts,max_edge_margin);
        }


        WaypointsProcessor(std::shared_ptr<valid_checker_t> valid_checker_ptr,
                           const std::vector<std::vector<double>> &pts,
                           double min_edge_margin = 1,double max_edge_margin = 8.0)
                           :WaypointsProcessor(valid_checker_ptr,DM(pts), min_edge_margin, max_edge_margin) {

        }
/*
        std::tuple<bool,std::vector<std::vector<double>>> optimize(std::shared_ptr<Path> path_ptr, int resolution=100,double ds_holder = 0.5){
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
            auto ones = DM::ones(1,4*steps+1);

            std::vector<size_t> index_vect(resolution+1);
            std::iota(index_vect.begin(),index_vect.end(),0);


            std::vector<double> inner_boundary_vec(resolution+1,2*max_edge_margin_);
            std::vector<double> outer_boundary_vec(resolution+1,-2*max_edge_margin_);

            std::for_each(std::execution::par,index_vect.begin(),index_vect.end(),[&](size_t idx){
                auto t_array = double(t(0,idx))*ones;
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
                    if(it>=8){
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
*/

        DM process(bool original_data = true) {

            DM points;
            std::vector<int> index_vec;
            if(original_data) {
                std::tie(points, index_vec) = split();
            }else{
                points = waypoints_;
                index_vec.resize(waypoints_.columns());
                std::iota(index_vec.begin(),index_vec.end(),0);
            }
            /*
            if(save_data)
                points_history.push_back(points);*/
            auto total_waypoints = points.columns();
            auto M = m_*DM::eye(total_waypoints);
            for(auto i:index_vec){
                M(i,i) = m_key_point_;
            }
            M(0,0) = 1e10;
            M(-1,-1)=1e10;

            for(auto idx = 0;idx<total_its;++idx){

                auto force_between_points = get_force_between_points(points);
                std::vector<std::vector<double>> pts;
                for(auto i=0;i<total_waypoints;++i){
                    pts.push_back({double(points(0,i)),double(points(1,i))});
                }

                std::vector<std::vector<double>> collision_force = valid_checker_ptr_->get_force(pts,max_edge_margin_);
                auto force_from_collision = DM(collision_force).T();
                auto total_force = force_between_points+force_from_collision;
                auto ax = DM::solve(M,total_force(0,Slice()).T());
                auto ay = DM::solve(M,total_force(1,Slice()).T());
                auto a = DM::horzcat({ax,ay});
                //std::cout<<"a=\n"<<a<<std::endl;
                //std::cout<<points.size();
                points = points + step_size_ * a.T();
                /*if(save_data)
                    points_history.push_back(points);*/
            }
            return points;
            /*
            DM v = points(Slice(),Slice(0,index_vec[1]));
            std::cout<<"Size: "<<points.rows()<<'\t'<<points.columns()<<std::endl;
            for(auto i =1;i<index_vec.size()-1;++i){
                v=DM::horzcat({v,points(Slice(),index_vec[i])});
            }
            v=DM::horzcat({v,points(Slice(),Slice(index_vec.back(),points.columns()))});
            std::cout<<"Size: "<<v.rows()<<'\t'<<v.columns()<<std::endl;
            return v;*/

        }

        std::vector<DM> get_history_data(){
            return points_history;
        }


    private:
        std::vector<DM> points_history;
        const double step_size_ = 1;
        const double waypoint_hold_distance = 30;
        const double k_distance_ =0.1;
        const double m_ = 10.0;
        const double m_key_point_ = 10.0;
        const int total_its = 200;
        //const int steps = 10;
        double max_edge_margin_{};
        double min_edge_margin_{};
        DM waypoints_;

        //DM dist_vec_;
        std::shared_ptr<valid_checker_t> valid_checker_ptr_;
        //std::shared_ptr<collision_field_t> collision_field_ptr_;
        std::shared_ptr<Path> path_ptr_;
        std::shared_ptr<Path> original_path_ptr_;


        std::tuple<DM,std::vector<int>> split() {
            DM dm;
            std::vector<int> index_vec;
            auto dm_diff = waypoints_(Slice(),Slice(1,waypoints_.columns())) - waypoints_(Slice(),Slice(0,-1));
            auto d = DM::sqrt(dm_diff(0,Slice())*dm_diff(0,Slice())+dm_diff(1,Slice())*dm_diff(1,Slice()));

            for(auto i=0;i<waypoints_.columns()-1;++i){
                index_vec.push_back(dm.columns());
                if(double(d(0,i))<=waypoint_hold_distance){
                    dm = DM::horzcat({dm,waypoints_(Slice(),i)});
                }else{
                    auto segment = std::ceil(double(d(0,i))/waypoint_hold_distance)+1;
                    auto new_dm = DM::linspace(waypoints_(Slice(),i).T(),waypoints_(Slice(),i+1).T(),segment);
                    //std::cout<<"new dm size: "<<new_dm.rows()<<'\t'<<new_dm.columns()<<std::endl;
                    dm = DM::horzcat({dm,new_dm.T()(Slice(),Slice(0,-1))});
                }

            }
            index_vec.push_back(dm.columns());

            dm = DM::horzcat({dm,waypoints_(Slice(),-1)});
            dm = DM::reshape(dm,2,dm.rows()*dm.columns()/2);


            return {dm,index_vec};
        }

        DM get_force_between_points(DM& points){
            DM force(2,points.columns());


            for(auto i=1;i<points.columns()-1;++i){
                force(Slice(),i) = k_distance_*(points(Slice(),i+1)-points(Slice(),i))+k_distance_*(points(Slice(),i-1)-points(Slice(),i));
            }
            force(Slice(),0) = k_distance_*(points(Slice(),1)-points(Slice(),0));
            force(Slice(),-1) = k_distance_*(points(Slice(),-2)-points(Slice(),-1));
            return force;
        }




    };
}


#endif //RACECARPLANNER_WAYPOINTS_PROCESSOR_HPP
