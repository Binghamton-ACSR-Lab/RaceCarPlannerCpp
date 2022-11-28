//
// Created by acsr on 11/10/22.
//

#ifndef CAR_CONTROL_CPP_GLOBAL_TRAJECTORY_HPP
#define CAR_CONTROL_CPP_GLOBAL_TRAJECTORY_HPP

#include "casadi/casadi.hpp"
#include "path.hpp"

using casadi::DM;
using casadi::Function;

namespace acsr {

    struct GlobalTrajectory {

        GlobalTrajectory() = default;

        GlobalTrajectory(std::shared_ptr<Path> path_ptr,const DM& dm_dt, const DM& dm_state, const DM& dm_control):path_ptr_(path_ptr) {
            assert(dm_dt.is_vector());
            DM control;
            if(dm_state.columns()==dm_control.columns()){
                control = dm_control;
            }else if(dm_state.columns()==dm_control.columns()+1){
                control = DM::zeros(dm_control.rows(),dm_control.columns()+1);
                control(Slice(),Slice(0,-1)) = dm_control(Slice(),Slice());

            }else{
                std::cerr<<"sizes of state and control not consistent"<<std::endl;
                return;
            }

            if(dm_dt->size() == dm_state.columns()){
                assert(double(dm_dt(0))<0.0001);
                dt_ = DM::cumsum(dm_dt);
            }else {
                assert(dm_dt->size()+1 == dm_state.columns());
                dt_ = DM::zeros(1, dm_dt->size() + 1);
                dt_(Slice(1, int(dt_->size()))) = DM::cumsum(dm_dt);
            }


            auto tau = dm_state(0,Slice());
            //auto tau = path_ptr->xy2t(dm_state(Slice(0,2),Slice()),true);
            auto tau_vec = std::vector<double>(tau->begin(),tau->end());
            auto dt_vec = std::vector<double>(dt_->begin(),dt_->end());


            f_tau_to_t = interpolant("tau_to_t","linear",std::vector<std::vector<double>>{tau_vec},dt_vec);
            //f_t_to_tau = interpolant("t_to_t","linear",std::vector<std::vector<double>>{dt_vec},tau_vec);
            f_t_to_state = interpolant("t_to_state","linear",std::vector<std::vector<double>>{dt_vec},std::vector<double>{dm_state->begin(),dm_state->end()});
            f_t_to_control = interpolant("t_to_control","linear",std::vector<std::vector<double>>{dt_vec},std::vector<double>{control->begin(),control->end()});
        }

        std::pair<DM,DM> get_reference_cartesian(const DM& start_state, int horizon,double dt,bool start_state_cartesian=true){
            DM tau0;
            if(start_state_cartesian) {
                tau0 = path_ptr_->xy2t(start_state(Slice(0, 2)));
            }else{
                tau0 = start_state(0);
            }
            auto t0 = f_tau_to_t(tau0)[0];
            auto t_vec = DM::linspace(t0,t0 + horizon*dt,horizon+1).T();
            auto state = f_t_to_state(t_vec)[0];
            auto xy = path_ptr_->f_tn_to_xy(DMVector{state(0,Slice()),state(1,Slice())})[0];
            state(Slice(0,2),Slice()) = xy;
            return std::make_pair(state,f_t_to_control(t_vec)[0]);
        }

        std::pair<DM,DM> get_reference_arc(const DM& start_state, int horizon,double dt,bool start_start_cartesian=true){
            DM tau0;
            if(start_start_cartesian) {
                tau0 = path_ptr_->xy2t(start_state(Slice(0, 2)));
            }else{
                tau0 = start_state(0);
            }
            auto t0 = f_tau_to_t(tau0)[0];
            auto t_vec = DM::linspace(t0,t0 + horizon*dt,horizon+1);
            return std::make_pair(f_t_to_state(t_vec)[0],f_t_to_control(t_vec)[0]);
        }



        std::shared_ptr<Path> get_path() const{
            return path_ptr_;
        }

    public:
        Function f_tau_to_t;
        Function f_t_to_tau;
        Function f_t_to_state;
        Function f_t_to_control;


    private:
        std::shared_ptr<Path> path_ptr_;
        DM dt_;
    };

}

#endif //CAR_CONTROL_CPP_GLOBAL_TRAJECTORY_HPP
