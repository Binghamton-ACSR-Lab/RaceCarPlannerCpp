//
// Created by acsr on 11/16/22.
//

#ifndef RACECARPLANNER_BICYCLE_KINEMATIC_HPP
#define RACECARPLANNER_BICYCLE_KINEMATIC_HPP

#include "utility.hpp"


namespace acsr {
    struct BicycleKinematic {
        BicycleKinematic(const param_t &param) {
            lf_ = param.at("lf");
            lr_ = param.at("lr");
        }


        template<class T, class op=std::function<T(const T &)>>
        T update(const T &x, const T &u, op f = nullptr) {
            auto phi = x(2, all);
            auto vx = x(3, all);

            auto delta = u(0, all);
            auto d = u(1, all);

            auto pt_x_dot = vx * T::cos(phi);
            auto pt_y_dot = vx * T::sin(phi);
            auto phi_dot = vx / (lf_ + lr_) * T::tan(delta);
            T v_dot;
            if (f == nullptr)
                v_dot = 7.9 * d;
            else
                v_dot = f(x);
            return T::vertcat({pt_x_dot, pt_y_dot, phi_dot, v_dot});
        }

        constexpr static unsigned int nx = 4;
        constexpr static unsigned int nu = 2;

    protected:
        const casadi::Slice all = casadi::Slice();
        double lf_;
        double lr_;

    };
}

#endif //RACECARPLANNER_BICYCLE_KINEMATIC_HPP
