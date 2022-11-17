//
// Created by acsr on 9/7/22.
//

#ifndef RACECARPLANNER_DYNAMICS_HPP
#define RACECARPLANNER_DYNAMICS_HPP

#include <casadi/casadi.hpp>
#include "track.hpp"

namespace acsr{
    using param_t = std::map<std::string,double>;



    template<class front_tire_model_t,class rear_tire_model_t>
    class BicycleDynamicsByParametricArc{
    public:
        BicycleDynamicsByParametricArc()=default;
        BicycleDynamicsByParametricArc(const param_t& param,std::shared_ptr<Track> track,std::shared_ptr<front_tire_model_t> front_tire_model,std::shared_ptr<rear_tire_model_t> rear_tire_model)
            :ptr_track{track},ptr_front_tire_model{front_tire_model},ptr_rear_tire_model{rear_tire_model}{
            lf = param.at("lf");
            lr = param.at("lr");
            m = param.at("m");
            Iz = param.at("Iz");

            cd = param.at("Cd");
            cm0 = param.at("Cm0");
            cm1 = param.at("Cm1");
            cm2 = param.at("Cm2");
            cbf = param.at("Cbf");
            cbr = param.at("Cbr");

        }

        template<class T>
        T updata(T& x,T& u){
            //x = [t,n,phi,vx,vy,omega]
            //u = [delta,d]

            auto N = x.columns()-1;
            auto t = x(0,Slice(0,N));
            auto n = x(1,Slice(0,N));
            auto phi = x(2,Slice(0,N));
            auto vx = x(3,Slice(0,N));
            auto vy = x(4,Slice(0,N));
            auto omega = x(5,Slice(0,N));
            auto delta = x(6,Slice(0,N));

            auto delta_dot = u(0,Slice());
            auto throttle = u(1,Slice());
            auto front_brake = u(2,Slice());
            auto rear_brake = u(3,Slice());

            auto kappa = ptr_track->f_kappa(t)[0];
            auto phi_c = ptr_track->f_phi(t)[0];
            auto tangent_vec = ptr_track->f_tangent_vec(t)[0];
            auto tangent_vec_norm = T::sqrt(tangent_vec(0,Slice())*tangent_vec(0,Slice())+tangent_vec(1,Slice())*tangent_vec(1,Slice()));
            auto Fz = casadi::DM::ones(2) * m / 2 * 9.81;
            //need work

            auto alphaf = -T::atan2(omega*lf + vy, vx) + delta;
            auto alphar = T::atan2(omega*lr - vy,vx);

            auto fx_f = cd * throttle - cm0 - cm1 * vx * vx - cm2 * vy * vy - cbf * front_brake;
            auto fx_r = cd * throttle - cm0 - cm1 * vx * vx - cm2 * vy * vy - cbr * rear_brake;

            auto fy_f = ptr_front_tire_model->get_lateral_force(alphaf, double(Fz(0)));
            auto fy_r = ptr_rear_tire_model->get_lateral_force(alphar, double(Fz(1)));

            auto t_dot = (vx*T::cos(phi-phi_c)-vy*T::sin(phi-phi_c))/(tangent_vec_norm*(1-n*kappa));
            auto n_dot = vx*T::sin(phi-phi_c)+vy*T::cos(phi-phi_c);
            auto phi_dot = omega;

            auto vx_dot = 1/m * (fx_r + fx_f*T::cos(delta) - fy_f*T::sin(delta) + m*vy*omega);
            auto vy_dot = 1/m * (fy_r + fx_f*T::sin(delta) + fy_f*T::cos(delta) - m*vx*omega);
            auto omega_dot = 1/Iz * (fy_f*lf*T::cos(delta) + fx_f*lf*T::sin(delta) - fy_r*lr);

            auto dot_x = T::vertcat({
                    t_dot,
                    n_dot,
                    phi_dot,
                    vx_dot,
                    vy_dot,
                    omega_dot,
                    delta_dot
            });
            return dot_x;
        }

        constexpr unsigned nx(){
            return 7;
        }

        constexpr unsigned nu(){
            return 4;
        }

    protected:
        double lr,lf;
        double m,Iz;
        double cd,cm0,cm1,cm2,cbf,cbr;



        std::shared_ptr<Track> ptr_track;
        std::shared_ptr<front_tire_model_t> ptr_front_tire_model;
        std::shared_ptr<rear_tire_model_t> ptr_rear_tire_model;
    };


}

#endif //RACECARPLANNER_DYNAMICS_HPP
