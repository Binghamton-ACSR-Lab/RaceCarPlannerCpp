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
    template<class valid_checker_t>
    class WaypointsProcessor {

    public:


        WaypointsProcessor(std::shared_ptr<valid_checker_t> valid_checker_ptr,
                           DM pts)
                           :valid_checker_ptr_(valid_checker_ptr){

            //format waypoints
            if (pts.rows() != 2) {
                assert(pts.columns() == 2);
                pts = pts.T();
            }
            waypoints_ = pts;

        }


        WaypointsProcessor(std::shared_ptr<valid_checker_t> valid_checker_ptr,
                           const std::vector<std::vector<double>> &pts)
                           :WaypointsProcessor(valid_checker_ptr,DM(pts)) {

        }

        DM process() {

            //construct new waypoints since the spacing of original waypoints may be large.
            DM original_points;
            std::vector<int> index_vec;
            int total_waypoints;
            bool success = false;
            int segment = 3;
            while (!success) {
                success=true;
                std::tie(original_points, index_vec) = split(segment++);
                total_waypoints = original_points.columns();
                for(auto i=1;i<total_waypoints-1;++i) {
                    if (!valid_checker_ptr_->valid(double(original_points(0,i)),double(original_points(1,i)))) {
                        success = false;
                        break;
                    }
                }
                if(segment>100){
                    std::cout<<"Original waypoints are invalid\n";
                    return DM{};
                }
            }





            //mass matrix, the mass of waypoint is different with others
            auto M = m_*DM::eye(total_waypoints-2);
            for(auto i=1;i<index_vec.size()-1;++i){
                M(index_vec[i]-1,index_vec[i]-1) = m_key_point_;
            }

            success = false;
            double step_size = 1.0;
            DM points;
            bool initial_iteration = true;
            while (!success) {
                points = original_points;
                for (auto idx = 0; idx < total_its; ++idx) {
                    success = true;
                    auto force_between_points = get_force_between_points(points);

                    std::vector<std::vector<double>> pts;
                    for (auto i = 1; i < total_waypoints - 1; ++i) {
                        pts.push_back({double(points(0, i)), double(points(1, i))});
                    }

                    std::vector<std::vector<double>> collision_force = valid_checker_ptr_->get_force(pts);
                    auto force_from_collision = DM(collision_force).T();
                    auto total_force = force_between_points(Slice(), Slice(1, -1)) + force_from_collision;
                    auto ax = DM::solve(M, total_force(0, Slice()).T());
                    auto ay = DM::solve(M, total_force(1, Slice()).T());
                    auto a = DM::horzcat({ax, ay});
                    points(Slice(), Slice(1, -1)) = points(Slice(), Slice(1, -1)) + step_size * a.T();


                    for(auto i=1;i<points.columns()-1;++i) {
                        if(!valid_checker_ptr_->valid(double(points(0,i)),double(points(1,i)))){
                            success = false;
                            step_size/=2;
                            break;
                        }
                    }
                    if(!success)break;
                }
            }
            return points;
        }



    private:
        const double k_hooke_ =1.0;
        const double m_ = 10.0;
        const double m_key_point_ = 20.0;
        const int total_its = 100;
        DM waypoints_;
        std::shared_ptr<valid_checker_t> valid_checker_ptr_;


        std::tuple<DM,std::vector<int>> split(int segment) {
            DM dm;
            std::vector<int> index_vec;
            auto dm_diff = waypoints_(Slice(),Slice(1,waypoints_.columns())) - waypoints_(Slice(),Slice(0,-1));
            auto d = DM::sqrt(dm_diff(0,Slice())*dm_diff(0,Slice())+dm_diff(1,Slice())*dm_diff(1,Slice()));

            auto d_sum = double(DM::sum2(d));
            auto hold_distance = d_sum/(segment*(waypoints_.columns()-1));

            for(auto i=0;i<waypoints_.columns()-1;++i){
                index_vec.push_back(dm.columns());
                if(double(d(0,i))<=hold_distance){
                    dm = DM::horzcat({dm,waypoints_(Slice(),i)});
                }else{
                    auto segment = std::ceil(double(d(0,i))/hold_distance)+1;
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
                force(Slice(),i) = k_hooke_*(points(Slice(),i+1)-points(Slice(),i))+k_hooke_*(points(Slice(),i-1)-points(Slice(),i));
            }
            force(Slice(),0) = k_hooke_*(points(Slice(),1)-points(Slice(),0));
            force(Slice(),-1) = k_hooke_*(points(Slice(),-2)-points(Slice(),-1));
            return force;
        }




    };
}


#endif //RACECARPLANNER_WAYPOINTS_PROCESSOR_HPP
