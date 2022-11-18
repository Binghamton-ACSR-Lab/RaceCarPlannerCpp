//
// Created by acsr on 11/16/22.
//

#ifndef RACECARPLANNER_GLOBAL_PLANNER_BICYCLE_KINEMATIC_ARC_HPP
#define RACECARPLANNER_GLOBAL_PLANNER_BICYCLE_KINEMATIC_ARC_HPP

#include "utility.hpp"
#include "path.hpp"
#include "global_planner.hpp"
#include "bicycle_kinematic_arc.hpp"

namespace acsr {


    class BicycleKinematicOptimizer {
    public:

        BicycleKinematicOptimizer(std::shared_ptr<Path> path_ptr,
                                           double path_width,
                                           const param_t &vehicle_params)
                : vehicle_params_(vehicle_params), path_width_(path_width),
                  path_ptr_(path_ptr),dynamics_{vehicle_params,path_ptr} {
            state_low_boundary_=std::vector<double>{0.0,0.0,0.0,0.1};
            state_up_boundary_=std::vector<double>{0.0,0.0,0.0,inf};
            control_low_boundary_=std::vector<double>{vehicle_params_.at("delta_min"),vehicle_params_.at("d_min")};
            control_up_boundary_=std::vector<double>{vehicle_params_.at("delta_max"),vehicle_params_.at("d_max")};


        }

        void make_plan(int horizon,double tau0, double s = -1,double n0 = 0,double phi0=-6, double v0 = 0.15) {
            auto goal=[this](const MX& dt,const MX& x,const MX& u){
                auto n_obj = MX::exp(100*(x(1,Slice())/(path_width_/2)-1)) + MX::exp(100*(x(1,Slice())/(-path_width_/2)-1));
                return 10*MX::sum2(dt) + MX::sum2(n_obj);
            };
            if(phi0<pi){
                phi0 = double(path_ptr_->f_phi(DM{tau0})[0]);
            }
            if(s<0.0){
                s = path_ptr_->get_max_length()-double(path_ptr_->t_to_s_lookup(DM{tau0})[0]);
            }
            auto s_vec = DM::linspace(0,s,horizon+1);
            auto t_vec = path_ptr_->s_to_t_lookup(s_vec)[0];
            t_vec = DM::reshape(t_vec,1,-1);

            auto x0 = DM::vertcat({tau0,n0,phi0,v0});
            auto return_value = global_planner_.make_plan(dynamics_,goal,horizon,x0,t_vec,state_low_boundary_,state_up_boundary_,control_low_boundary_,control_up_boundary_);

        }




    private:
        //constexpr static int nx = 4;
        //constexpr static int nu = 2;
        BicycleKinematicArc dynamics_;
        AcsrGlobalPlanner<BicycleKinematicArc> global_planner_;

        std::shared_ptr<Path> path_ptr_;
        double path_width_{};
        param_t vehicle_params_, track_params_;

        std::vector<double> state_low_boundary_,state_up_boundary_,control_low_boundary_,control_up_boundary_;

        //std::shared_ptr<valid_checker_t> valid_checker_;



    };

}


#endif //RACECARPLANNER_GLOBAL_PLANNER_BICYCLE_KINEMATIC_ARC_HPP
