//
// Created by acsr on 11/9/22.
//

#ifndef CAR_CONTROL_CPP_REMOTE_CONTROLLER_HPP_
#define CAR_CONTROL_CPP_REMOTE_CONTROLLER_HPP_

#include "path.hpp"
#include "../unity_controller/include/websocket_server.hpp"
#include <casadi/casadi.hpp>
#include "global_trajectory.hpp"

using json = nlohmann::json;
using std::string;
using casadi::DM;

namespace acsr {

    struct KinematicModelController {

        KinematicModelController() = default;

        KinematicModelController(std::shared_ptr<GlobalTrajectory> gloabal_path_ptr, int horizon, double dt):global_path_ptr_(gloabal_path_ptr),horizon_(horizon),dt_(dt){
            x = opti.variable(nx,horizon_+1);
            u = opti.variable(nu,horizon_);

        }

        std::pair<DM,DM> make_plan(const DM& x0){
            DM ref_state,ref_control;
            auto x_dot = dynamics_model_cartesian<MX>(x(Slice(),Slice(0,-1)),u);
            std::tie(ref_state,ref_control) = global_path_ptr_->get_reference(x0,horizon_,dt_);


            opti.minimize(5 * MX::dot((x(0, all) - ref_state(0,all)),(x(0, all) - ref_state(0,all)))
                          + 5 * MX::dot((x(1, all) - ref_state(1,all)),(x(1, all) - ref_state(1,all))));

            //dynamics
            opti.subject_to(x(Slice(),Slice(1,horizon_+1))==x(Slice(),Slice(0,-1))+dt_*x_dot);

            //initial value
            opti.subject_to(x(Slice(),0)==x0);

            //state boundary
            opti.subject_to(opti.bounded(v_min_,x(3,all),v_max_));

            //control boundary
            opti.subject_to(opti.bounded(d_min_,u(1,all),d_max_));
            opti.subject_to(opti.bounded(delta_min_,u(0,all),delta_max_));

            //initial guess
            opti.set_initial(x, ref_state);
            opti.set_initial(u, ref_control(all,Slice(0,-1)));

            try{
                auto sol = opti.solve();
                auto sol_x = sol.value(x);
                auto sol_u = sol.value(u);
                return std::make_pair(sol_x,sol_u);

            }catch(CasadiException& e) {
                std::cout<<e.what()<<std::endl;
                std::cout<<"Solve Optimal Problem Fails\n";
                return std::make_pair(DM{},DM{});
            }
        }

        json operator()(const json &value) {
            double phi0 = value.at("psi");
            double v0 = value.at("speed");
            double pt_x0 = value.at("x");
            double pt_y0 = value.at("y");

            auto x0 = DM::vertcat({pt_x0,pt_y0,phi0,v0});
            DM sol_x,sol_u;
            std::tie(sol_x,sol_u) = make_plan(x0);
            if(sol_x->size()==0){
                return R"({"foo": "bar"})"_json;
            }else{

            }


        }

    public:
        static constexpr int nx = 4;
        static constexpr int nu = 2;

    private:
        const Slice all = casadi::Slice();
        DM dm_waypoints;
        std::shared_ptr<GlobalTrajectory> global_path_ptr_;
        //std::shared_ptr<Path> path_ptr_;
        double d_min_, d_max_, delta_min_, delta_max_, lf_, lr_, mass_, Iz,v_min_,v_max_;


        casadi::Opti opti;
        MX x;
        MX u;
        double dt_;
        int horizon_;

        template<class T,class op=std::function<T(const T&)>>
        T dynamics_model_arc(const T& x, const T& u,op f =nullptr){
            auto t = x(0,all);
            auto n = x(1,all);
            auto phi = x(2,all);
            auto vx = x(3,all);

            auto delta = u(0,all);
            auto d = u(1,all);

            auto path_ptr = global_path_ptr_->get_path();

            auto kappa = path_ptr->f_kappa(t);
            auto phi_c = path_ptr->f_phi(t);
            auto tangent_vec = path_ptr->f_tangent_vec(t);

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

        template<class T,class op=std::function<T(const T&)>>
        T dynamics_model_cartesian(const T& x, const T& u,op f =nullptr){
            auto phi = x(2,all);
            auto vx = x(3,all);

            auto delta = u(0,all);
            auto d = u(1,all);

            auto pt_x_dot = vx*T::cos(phi);
            auto pt_y_dot = vx*T::sin(phi);
            auto phi_dot = vx/(lf_+lr_) * T::tan(delta);
            T v_dot;
            if(f==nullptr)
                v_dot = 7.9*d;
            else
                v_dot = f(x);
            return T::vertcat({pt_x_dot,pt_y_dot,phi_dot,v_dot});
        }

    };
}

#endif //CAR_CONTROL_CPP_REMOTE_CONTROLLER_HPP_
