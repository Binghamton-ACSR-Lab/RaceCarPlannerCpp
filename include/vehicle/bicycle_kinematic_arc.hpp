//
// Created by acsr on 11/16/22.
//

#ifndef RACECARPLANNER_BICYCLE_KINEMATIC_ARC_HPP
#define RACECARPLANNER_BICYCLE_KINEMATIC_ARC_HPP
#include "utility.hpp"
#include "path.hpp"
namespace acsr {


    struct BicycleKinematicArc {

        BicycleKinematicArc(const param_t &param, std::shared_ptr<Path> path_ptr):path_ptr_(path_ptr){
            if(param.find("wheel_base")!=param.end())
                wheel_base_ = param.at("wheel_base");
            else
                wheel_base_ = param.at("lf")+param.at("lr");
        }


        MX update(const MX& x, const MX& u){
            //auto t = x(0,all);
            //auto n = x(1,all);
            //auto phi = x(2,all);
            //auto vx = x(3,all);

            //auto delta = u(0,all);
            //auto d = u(1,all);

            //auto path_ptr = global_path_ptr_->get_path();

            auto kappa = path_ptr_->f_kappa(x(0))[0];
            auto phi_c = path_ptr_->f_phi(x(0))[0];
            auto tangent_vec = path_ptr_->f_tangent_vec(x(0))[0];

            auto t_dot = x(3)*MX::cos(x(2)-phi_c+u(0))/(MX::norm_2(tangent_vec)*(1.0-x(1)*kappa));
            auto n_dot = x(3)*MX::sin(x(2)-phi_c+u(0));
            auto phi_dot = x(3)/wheel_base_ * MX::tan(u(0));
            //auto v_dot = ;
            return MX::vertcat({t_dot,n_dot,phi_dot,7.9*u(1)});
        }

        constexpr static unsigned int nx = 4;
        constexpr static unsigned int nu = 2;

    private:
        const casadi::Slice all = casadi::Slice();
        double wheel_base_{};

        std::shared_ptr<Path> path_ptr_;

    };
}

#endif //RACECARPLANNER_BICYCLE_KINEMATIC_ARC_HPP
