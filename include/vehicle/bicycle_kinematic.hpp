//
// Created by acsr on 11/16/22.
//

#ifndef RACECARPLANNER_BICYCLE_KINEMATIC_HPP
#define RACECARPLANNER_BICYCLE_KINEMATIC_HPP

#include "utility.hpp"


namespace acsr {


    struct BicycleKinematic {
        BicycleKinematic(const param_t &param) {
            if(param.find("wheel_base")!=param.end())
                wheel_base_ = param.at("wheel_base");
            else
                wheel_base_ = param.at("lf")+param.at("lr");
        }


        template<class T>
        T operator()(const T &x, const T &u) {
            auto phi = x(2, all);
            auto vx = x(3, all);

            auto delta = u(0, all);
            auto d = u(1, all);

            auto pt_x_dot = vx * T::cos(phi);
            auto pt_y_dot = vx * T::sin(phi);
            auto phi_dot = vx / wheel_base_ * T::tan(delta);
            T v_dot = 7.9 * d;
            return T::vertcat({pt_x_dot, pt_y_dot, phi_dot, v_dot});
        }

        template<class T,class OP>
        T operator()(const T &x, const T &u,OP f) {
            auto phi = x(2, all);
            auto vx = x(3, all);

            auto delta = u(0, all);
            auto d = u(1, all);

            auto pt_x_dot = vx * T::cos(phi);
            auto pt_y_dot = vx * T::sin(phi);
            auto phi_dot = vx / wheel_base_ * T::tan(delta);
            T v_dot = f(x,u);
            return T::vertcat({pt_x_dot, pt_y_dot, phi_dot, v_dot});
        }

        constexpr static unsigned int nx = 4;
        constexpr static unsigned int nu = 2;

    protected:
        const casadi::Slice all = casadi::Slice();
        //double lf_;
        //double lr_;
        double wheel_base_{};

    };
}

#endif //RACECARPLANNER_BICYCLE_KINEMATIC_HPP
