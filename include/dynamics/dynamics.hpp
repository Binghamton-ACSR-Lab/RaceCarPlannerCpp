//
// Created by xli185 on 12/16/22.
//

#ifndef RACECARPLANNER_DYNAMICS_HPP
#define RACECARPLANNER_DYNAMICS_HPP

#include <functional>
#include <casadi/casadi.hpp>
#include <nlohmann/json.hpp>
#include "path.hpp"
#include "utility.hpp"
#include <algorithm>

using casadi::Slice;
using nlohmann::json;

namespace acsr {

    template<int nx, int nu>
    struct Dynamics {
    public:
        static constexpr int nx_ = nx;
        static constexpr int nu_ = nu;

        using StateType = std::array<double,nx>;
        using ControlType = std::array<double,nu>;
        //using CellType = std::array<int,CELL_DIMENSION>;
        //using TreeIndexType = std::array<double,TREE_INDEX_DIMENSION>;

        Dynamics() = default;

        virtual double euclidean_distance(const StateType &x1, const StateType &x2){
            return sqrt((x1[0]-x2[0])*(x1[0]-x2[0])+(x1[1]-x2[1])*(x1[1]-x2[1]));
        }
    };


    struct BicycleKinematic : Dynamics<4, 2> {

//        using StateType = std::array<double,4>;
//        using ControlType = std::array<double,2>;
//        using CellType = std::array<int,CELL_DIMENSION>;
//        using TreeIndexType = std::array<double,TREE_INDEX_DIMENSION>;

        BicycleKinematic() = default;

        BicycleKinematic(const json &params) {
            if(params.contains("wheel_base"))
                wheel_base_ = params.at("wheel_base");
            else
                wheel_base_ = double(params.at("lf"))+double(params.at("lr"));

            cd_=params.contains("cd")?double(params.at("cd")):7.9;
            cm0_=params.contains("cm0")?double(params.at("cm0")):0.0;
            cm1_=params.contains("cm1")?double(params.at("cm1")):0.0;
            cm2_=params.contains("cm2")?double(params.at("cm2")):0.0;

            v_min_=params.at("v_min");
            v_max_=params.at("v_max");
            a_min_=params.at("a_min");
            a_max_=params.at("a_max");
            delta_min_=params.at("delta_min");
            delta_max_=params.at("delta_max");
        }

        template<class T>
        T model_arc(const T &x, const T &u, std::shared_ptr<Path> path_ptr) {
            auto t = x(0, all);
            auto n = x(1, all);
            auto phi = x(2, all);
            auto vx = x(3, all);

            auto delta = u(0, all);
            auto d = u(1, all);

            auto kappa = path_ptr->f_kappa(t);
            auto phi_c = path_ptr->f_phi(t);
            auto tangent_vec = path_ptr->f_tangent_vec(t);

            auto t_dot = vx * T::cos(phi - phi_c + delta) / (T::norm_2(tangent_vec) * (1.0 - n * kappa));
            auto n_dot = vx * T::sin(phi - phi_c + delta);
            auto phi_dot = vx / wheel_base_ * T::tan(delta);
            T v_dot = cd_ * d - cm0_ - cm1_ * vx - cm2_ * vx * vx;
            return T::vertcat({t_dot, n_dot, phi_dot, v_dot});
        }

        template<class T>
        T model(const T &x, const T &u) {
            auto phi = x(2, all);
            auto vx = x(3, all);

            auto delta = u(0, all);
            auto d = u(1, all);

            auto pt_x_dot = vx * T::cos(phi);
            auto pt_y_dot = vx * T::sin(phi);
            auto phi_dot = vx / wheel_base_ * T::tan(delta);
            T v_dot = cd_ * d - cm0_ - cm1_ * vx - cm2_ * vx * vx;

            return T::vertcat({pt_x_dot, pt_y_dot, phi_dot, v_dot});
        }

        template<class MAP>
        std::pair<bool,StateType > forward(std::shared_ptr<MAP> map, const StateType &init_state, const ControlType &control, double dt,int steps){
            auto result_state = init_state;
            for(auto i=0;i<steps;++i){
                result_state[0] = result_state[0]+result_state[3]*cos(control[0]+result_state[2])*dt;
                result_state[1] = result_state[1]+result_state[3]*sin(control[0]+result_state[2])*dt;
                result_state[2] = result_state[2]+result_state[3]*tan(control[0])/wheel_base_*dt;
                auto a = cd_ * control[1] - cm0_ - cm1_ * result_state[3] - cm2_ * result_state[3] * result_state[3];
                result_state[3] = std::clamp(result_state[3]+a*dt,v_min_,v_max_);

                if(!map->validate(result_state)){
                    return {false,result_state};
                }
            }
            return {true,result_state};
        }

        template<class MAP>
        std::pair<bool,StateType > backward(std::shared_ptr<MAP> map, const StateType &init_state, const ControlType &control, double dt,int steps){
            auto result_state = init_state;
            for(auto i=0;i<steps;++i){
                result_state[0] = result_state[0]-result_state[3]*cos(control[0]+result_state[2])*dt;
                result_state[1] = result_state[1]-result_state[3]*sin(control[0]+result_state[2])*dt;
                result_state[2] = result_state[2]-result_state[3]*tan(control[0])/wheel_base_*dt;
                result_state[3] = result_state[3]+control[1]*dt;
                auto a = cd_ * control[1] - cm0_ - cm1_ * result_state[3] - cm2_ * result_state[3] * result_state[3];
                result_state[3] = std::clamp(result_state[3]+a*dt,v_min_,v_max_);

                if(!map->validate(result_state)){
                    return {false,result_state};
                }
            }
            return {true,result_state};
        }

        template<class MAP>
        StateType random_state(std::shared_ptr<MAP> map){
            StateType state;
            auto& rnd = acsr::Random::get_instance();
            auto low_bounday = map->get_low_boundary();
            auto up_bounday=map->get_up_boundary();
            state[0] = rnd.random_value(low_bounday[0],up_bounday[0]);
            state[1] = rnd.random_value(low_bounday[1],up_bounday[1]);
            state[2] = rnd.random_value(-M_PI,M_PI);
            state[3] = rnd.random_value(v_min_,v_max_);
            return state;
        }

        ControlType random_control(){
            auto& rnd = acsr::Random::get_instance();
            return ControlType{rnd.random_value(delta_min_,delta_max_),rnd.random_value(a_min_,a_max_)};
        }



    private:
        const Slice all = casadi::Slice();
        double cd_, cm0_, cm1_, cm2_;
        double wheel_base_;
        double delta_min_,delta_max_;
        double v_min_,v_max_;
        double a_min_,a_max_;

    };

    struct BicycleKinematicNoAcc : Dynamics<3,2>{
        BicycleKinematicNoAcc()=default;

        BicycleKinematicNoAcc(const json &params) {
            if(params.contains("wheel_base"))
                wheel_base_ = params.at("wheel_base");
            else
                wheel_base_ = double(params.at("lf"))+double(params.at("lr"));

            v_min_=params.at("v_min");
            v_max_=params.at("v_max");

            delta_min_=params.at("delta_min");
            delta_max_=params.at("delta_max");

        }

        template<class T>
        T model_arc(const T& x, const T& u, std::shared_ptr<Path> path_ptr){
            auto t = x(0,all);
            auto n = x(1,all);
            auto phi = x(2,all);

            auto delta = u(0,all);
            auto v = u(1,all);

            auto kappa = path_ptr->f_kappa(t);
            auto phi_c = path_ptr->f_phi(t);
            auto tangent_vec = path_ptr->f_tangent_vec(t);

            auto t_dot = v*T::cos(phi-phi_c+delta)/(T::norm_2(tangent_vec)*(1.0-n*kappa));
            auto n_dot = v*T::sin(phi-phi_c+delta);
            auto phi_dot = v/wheel_base_ * T::tan(delta);

            return T::vertcat({t_dot,n_dot,phi_dot});
        }

        template<class T>
        T model(const T& x, const T& u){
            auto phi = x(2,all);

            auto delta = u(0,all);
            auto v = u(1,all);

            auto pt_x_dot = v*T::cos(phi);
            auto pt_y_dot = v*T::sin(phi);
            auto phi_dot = v/wheel_base_ * T::tan(delta);

            return T::vertcat({pt_x_dot,pt_y_dot,phi_dot});
        }

        template<class MAP>
        std::pair<bool,StateType > forward(std::shared_ptr<MAP> map, const StateType &init_state, const ControlType &control, double dt,int steps){
            auto result_state = init_state;
            for(auto i=0;i<steps;++i){
                result_state[0] = result_state[0]+control[1]*cos(control[0]+result_state[2])*dt;
                result_state[1] = result_state[1]+control[1]*sin(control[0]+result_state[2])*dt;
                result_state[2] = result_state[2]+control[1]*tan(control[0])/wheel_base_*dt;

                if(!map->validate(result_state)){
                    return {false,result_state};
                }
            }
            return {true,result_state};
        }

        template<class MAP>
        std::pair<bool,StateType > backward(std::shared_ptr<MAP> map, const StateType &init_state, const ControlType &control, double dt,int steps){
            auto result_state = init_state;
            for(auto i=0;i<steps;++i){
                result_state[0] = result_state[0]-control[1]*cos(control[0]+result_state[2])*dt;
                result_state[1] = result_state[1]-control[1]*sin(control[0]+result_state[2])*dt;
                result_state[2] = result_state[2]-control[1]*tan(control[0])/wheel_base_*dt;

                if(!map->validate(result_state)){
                    return {false,result_state};
                }
            }
            return {true,result_state};
        }

        template<class MAP>
        StateType random_state(std::shared_ptr<MAP> map){
            StateType state;
            auto& rnd = acsr::Random::get_instance();
            auto low_bounday = map->get_low_boundary();
            auto up_bounday=map->get_up_boundary();
            state[0] = rnd.random_value(low_bounday[0],up_bounday[0]);
            state[1] = rnd.random_value(low_bounday[1],up_bounday[1]);
            state[2] = rnd.random_value(-M_PI,M_PI);
            return state;
        }

        ControlType random_control(){
            auto& rnd = acsr::Random::get_instance();
            return ControlType{rnd.random_value(delta_min_,delta_max_),rnd.random_value(v_min_,v_max_)};
        }

    private:
        double wheel_base_;
        double delta_min_,delta_max_;
        double v_min_,v_max_;
    };
}

#endif //RACECARPLANNER_DYNAMICS_H
