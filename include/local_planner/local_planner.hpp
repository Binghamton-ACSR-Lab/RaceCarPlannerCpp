//
// Created by acsr on 11/9/22.
//

#ifndef CAR_CONTROL_CPP_REMOTE_CONTROLLER_HPP_
#define CAR_CONTROL_CPP_REMOTE_CONTROLLER_HPP_

#include "path.hpp"
#include "include/websocket_server.hpp"
#include <casadi/casadi.hpp>
#include "global_trajectory.hpp"
#include "db_manager.hpp"

using json = nlohmann::json;
using std::string;
using casadi::DM;

namespace acsr {

    template<int nx,int nu>
    class AcsrLocalPlanner{
    public:
        static constexpr int nx_ = nx;
        static constexpr int nu_ = nu;
        AcsrLocalPlanner()=default;

        virtual std::pair<bool,std::pair<DM,DM>> make_plan(const DM& x0,bool print) =0;


        virtual std::string save(DbManager &db_manager, DM& x,DM& u,const std::vector<std::string>& x_headers={},const std::vector<std::string>& u_headers={}){
            //std::cout << "Save data to SQLite database file '" << db.getFilename().c_str() << "\n";
            const auto now = std::chrono::system_clock::now();
            time_t rawtime;
            struct tm * timeinfo;
            char buffer[40];
            time (&rawtime);
            timeinfo = localtime(&rawtime);
            strftime(buffer,sizeof(buffer),"_%d_%H_%M_%S",timeinfo);
            auto datatable_name = std::string(buffer);

            auto local_x_header = x_headers;
            auto local_u_header = u_headers;
            if(local_x_header.size()<nx_){
                for (int i = 0; i < nx_; ++i) {
                    local_x_header.push_back("x_"+std::to_string(i));
                }
            }
            if(local_u_header.size()<nu_){
                for (int i = 0; i < nu_; ++i) {
                    local_u_header.push_back("u_"+std::to_string(i));
                }
            }

            db_manager.add_sql("DROP TABLE IF EXISTS " + datatable_name);

            std::ostringstream os;
            os<<"create table if not exists "<<datatable_name<<"(id INTEGER PRIMARY KEY AUTOINCREMENT, dt REAL ";
            for (int i = 0; i < nx_; ++i) {
                os<<","<<local_x_header[i]<<" REAL";
            }

            for (int i = 0; i < nu_; ++i) {
                os<<","<<local_u_header[i]<<" REAL";
            }
            os<<")";
            db_manager.add_sql(os.str());

            std::ostringstream os1;
            os1 << "INSERT INTO "<<datatable_name<<" (dt";
            for (int i = 0; i < nx_; ++i) {
                os1<<","<<local_x_header[i];
            }
            for (int i = 0; i < nu_; ++i) {
                os1<<","<<local_u_header[i];
            }
            os1 <<") VALUES(";
            auto init_string = os1.str();
            //insert data rows
            auto N = x.columns()-1;
            int divider = N/10;
            std::vector<std::string> vec;
            for(auto i=0;i<N;++i){
                if(i%divider==0){
                    std::cout<<"Saving data... "<<i*10/divider<<"% finished"<<std::endl;
                }
                std::ostringstream os;
                os<<double(x(0,i));
                for(auto j=1;j<nx_;++j){
                    os<<","<<double(x(j,i));
                }
                for(auto j=0;j<nu_;++j){
                    os<<","<<double(u(0,i));
                }
                os<<")";
                vec.push_back(init_string+os.str());
            }

            std::ostringstream os2;
            os2 << "INSERT INTO "<<datatable_name<<" ("<<local_x_header[0];
            for (int i = 1; i < nx_; ++i) {
                os2<<","<<local_x_header[i];
            }
            os2 <<") VALUES("<<double(x(0,N));

            for(auto j=1;j<nx_;++j){
                os2<<","<<double(x(j,N));
            }
            os2<<")";
            vec.push_back(os2.str());
            db_manager.add_sql(vec.begin(),vec.end());

            std::cout<<"save to database finished\n";
            return datatable_name;
        }

    };


    class BicycleKinematicController : public AcsrLocalPlanner<4,2> {

    public:

        BicycleKinematicController() = default;

        BicycleKinematicController(std::shared_ptr<GlobalTrajectory> gloabal_path_ptr, const json& params,int horizon, double dt):global_path_ptr_(gloabal_path_ptr),horizon_(horizon),dt_(dt){

            auto& constraint = params.at("constraint");

            delta_min_ = constraint.at("delta_min");
            delta_max_ = constraint.at("delta_max");
            v_min_ = constraint.at("v_min");
            v_max_ = constraint.at("v_max");
            d_min_ = constraint.at("d_min");
            d_max_ = constraint.at("d_max");

            if(params.contains("wheel_base"))
                wheel_base_ = params.at("wheel_base");
            else
                wheel_base_ = double(params.at("lf"))+double(params.at("lr"));

            option_["max_iter"] = 3000;
            option_["tol"]=1e-9;
            option_["linear_solver"]="ma57";

            weight_ = DM::linspace(0,1,horizon_+1).T();
            //weight_ = DM::exp(DM::linspace(0,5,horizon_+1).T());
        }

        std::pair<bool,std::pair<DM,DM>> make_plan(const DM& x0, bool print = false){
            casadi::Opti opti;
            auto x = opti.variable(nx_,horizon_+1);
            auto u = opti.variable(nu_,horizon_);
            //MX x;
            //MX u;

            DM ref_state,ref_control;
            auto x_dot = dynamics_model_cartesian<MX>(x(Slice(),Slice(0,-1)),u);


            std::tie(ref_state,ref_control) = global_path_ptr_->get_reference_cartesian(x0,horizon_,dt_);

            //std::cout<<"x0:"<<x0<<std::endl;
            //std::cout<<"ref_x:\n"<<ref_state<<std::endl;
            //std::cout<<"ref_u:\n"<<ref_control<<std::endl;
            auto diff_x = weight_*(x(0, all) - ref_state(0,all));
            auto diff_y = weight_*(x(1, all) - ref_state(1,all));
            auto diff_phi = weight_*(x(2, all) - ref_state(2,all));
            auto diff_v = weight_*(x(3, all) - ref_state(3,all));
            auto steer_diff = u(0,Slice(1,horizon_))-u(0,Slice(0,-1));

            opti.minimize(0.5 * MX::dot(diff_x,diff_x)
                          + 0.5 * MX::dot(diff_y,diff_y)
                            + 0.5 * MX::dot(diff_phi,diff_phi)
                            + 5 * MX::dot(diff_v,diff_v)
                            + 0.1 * MX::dot(steer_diff,steer_diff));

            //dynamics
            opti.subject_to(x(Slice(),Slice(1,horizon_+1))==x(Slice(),Slice(0,-1))+dt_*x_dot);

            //initial value
            opti.subject_to(x(Slice(),0)==x0);

            //state boundary
            //opti.subject_to(opti.bounded(v_min_,x(3,all),v_max_));

            //control boundary
            opti.subject_to(opti.bounded(d_min_,u(1,all),d_max_));
            opti.subject_to(opti.bounded(delta_min_,u(0,all),delta_max_));

            //initial guess
            //opti.set_initial(x, ref_state);
            //opti.set_initial(u, ref_control(all,Slice(0,-1)));

            Dict casadi_option;
            casadi_option["print_time"]=print;
            if(print){
                option_["print_level"]=5;
            }else{
                option_["print_level"]=1;
            }

            opti.solver("ipopt", casadi_option, option_);


            try{
                auto sol = opti.solve();
                auto sol_x = sol.value(x);
                auto sol_u = sol.value(u);
                //std::cout<<"x:\n"<<sol_x<<std::endl;
                //std::cout<<"u:\n"<<sol_x<<std::endl;
                return std::make_pair(true,std::make_pair(sol_x,sol_u));

            }catch(CasadiException& e) {
                std::cout<<e.what()<<std::endl;
                std::cout<<"Solve Optimal Problem Fails\n";
                return std::make_pair(false,std::make_pair(DM{},DM{}));
            }
        }

    private:
        const Slice all = casadi::Slice();
        DM dm_waypoints;
        std::shared_ptr<GlobalTrajectory> global_path_ptr_;
        double d_min_, d_max_, delta_min_, delta_max_, wheel_base_, v_min_,v_max_;

        Dict option_;

        double dt_;
        int horizon_;

        DM weight_;


    public:
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
            auto phi_dot = vx/wheel_base_ * T::tan(delta);
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
            auto phi_dot = vx/wheel_base_ * T::tan(delta);
            T v_dot;
            if(f==nullptr)
                v_dot = 7.9*d;
            else
                v_dot = f(x);
            return T::vertcat({pt_x_dot,pt_y_dot,phi_dot,v_dot});
        }

    };

    class BicycleKinematicNoAccController : public AcsrLocalPlanner<3,2> {

    public:

        BicycleKinematicNoAccController() = default;

        BicycleKinematicNoAccController(std::shared_ptr<GlobalTrajectory> gloabal_path_ptr, const json& params,int horizon, double dt):
                global_path_ptr_(gloabal_path_ptr),horizon_(horizon),dt_(dt){
            x = opti.variable(nx_,horizon_+1);
            u = opti.variable(nu_,horizon_);
            auto& constraint = params.at("constraint");

            delta_min_ = constraint.at("delta_min");
            delta_max_ = constraint.at("delta_max");
            v_min_ = constraint.at("v_min");
            v_max_ = constraint.at("v_max");
//            d_min_ = constraint.at("d_min");
//            d_max_ = constraint.at("d_max");

            if(params.contains("wheel_base"))
                wheel_base_ = params.at("wheel_base");
            else
                wheel_base_ = double(params.at("lf"))+double(params.at("lr"));

        }

        std::pair<DM,DM> make_plan(const DM& x0){
            DM ref_state,ref_control;
            auto x_dot = dynamics_model_cartesian<MX>(x(Slice(),Slice(0,-1)),u);
            std::tie(ref_state,ref_control) = global_path_ptr_->get_reference_cartesian(x0,horizon_,dt_,true);


            opti.minimize(5 * MX::dot((x(0, all) - ref_state(0,all)),(x(0, all) - ref_state(0,all)))
                          + 5 * MX::dot((x(1, all) - ref_state(1,all)),(x(1, all) - ref_state(1,all))));

            //dynamics
            opti.subject_to(x(Slice(),Slice(1,horizon_+1))==x(Slice(),Slice(0,-1))+dt_*x_dot);

            //initial value
            opti.subject_to(x(Slice(),0)==x0);

            //state boundary


            //control boundary
            //opti.subject_to(opti.bounded(d_min_,u(1,all),d_max_));
            opti.subject_to(opti.bounded(v_min_,u(1,all),v_max_));
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

    private:
        const Slice all = casadi::Slice();
        DM dm_waypoints;
        std::shared_ptr<GlobalTrajectory> global_path_ptr_;
        double d_min_, d_max_, delta_min_, delta_max_, wheel_base_, v_min_,v_max_;


        casadi::Opti opti;
        MX x;
        MX u;
        double dt_;
        int horizon_;

        template<class T>
        T dynamics_model_arc(const T& x, const T& u){
            auto t = x(0,all);
            auto n = x(1,all);
            auto phi = x(2,all);
            //auto vx = x(3,all);

            auto delta = u(0,all);
            auto v = u(1,all);

            auto path_ptr = global_path_ptr_->get_path();

            auto kappa = path_ptr->f_kappa(t);
            auto phi_c = path_ptr->f_phi(t);
            auto tangent_vec = path_ptr->f_tangent_vec(t);

            auto t_dot = v*T::cos(phi-phi_c+delta)/(T::norm_2(tangent_vec)*(1.0-n*kappa));
            auto n_dot = v*T::sin(phi-phi_c+delta);
            auto phi_dot = v/wheel_base_ * T::tan(delta);

            return T::vertcat({t_dot,n_dot,phi_dot});
        }

        template<class T>
        T dynamics_model_cartesian(const T& x, const T& u){
            auto phi = x(2,all);
            //auto vx = x(3,all);

            auto delta = u(0,all);
            auto v = u(1,all);

            auto pt_x_dot = v*T::cos(phi);
            auto pt_y_dot = v*T::sin(phi);
            auto phi_dot = v/wheel_base_ * T::tan(delta);

            return T::vertcat({pt_x_dot,pt_y_dot,phi_dot});
        }

    };
}

#endif //CAR_CONTROL_CPP_REMOTE_CONTROLLER_HPP_
