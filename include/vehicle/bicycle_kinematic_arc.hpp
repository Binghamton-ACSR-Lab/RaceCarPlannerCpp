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
            lf_ = param.at("lf");
            lr_ = param.at("lr");
        }

        template<class T,class op=std::function<T(const T&)>>
        T update(const T& x, const T& u,op f =nullptr){
            auto t = x(0,all);
            auto n = x(1,all);
            auto phi = x(2,all);
            auto vx = x(3,all);

            auto delta = u(0,all);
            auto d = u(1,all);

            //auto path_ptr = global_path_ptr_->get_path();

            auto kappa = path_ptr_->f_kappa(t);
            auto phi_c = path_ptr_->f_phi(t);
            auto tangent_vec = path_ptr_->f_tangent_vec(t);

            auto t_dot = vx*T::cos(phi-phi_c+delta)/(T::norm_2(tangent_vec)*(1.0-n*kappa));
            auto n_dot = vx*T::sin(phi-phi_c+delta);
            auto phi_dot = vx/(lf_+lr_) * T::tan(delta);
            T v_dot;
            if(f==nullptr)
                v_dot = 7.9*d;
            else
                v_dot = f(x);
            return T::vertcat({t_dot,n_dot,phi_dot,v_dot});
        }

        constexpr static unsigned int nx = 4;
        constexpr static unsigned int nu = 2;

    private:
        const casadi::Slice all = casadi::Slice();
        double lf_;
        double lr_;
        std::shared_ptr<Path> path_ptr_;

    };
}

#endif //RACECARPLANNER_BICYCLE_KINEMATIC_ARC_HPP
