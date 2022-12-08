//
// Created by acsr on 9/27/22.
//

#ifndef RACECARPLANNER_GLOBAL_PLANNER_OPTIMIZER_HPP
#define RACECARPLANNER_GLOBAL_PLANNER_OPTIMIZER_HPP

#include "path.hpp"
#include <filesystem>
#include "dynamics.hpp"
#include "tire_model.hpp"
#include <SQLiteCpp/SQLiteCpp.h>
#include <boost/format.hpp>
#include <filesystem>
#include <execution>
#include "nlohmann/json.hpp"
#include "db_manager.hpp"

namespace acsr{
    using param_t = std::map<std::string,double>;
    using json = nlohmann::json;


    template<int nx,int nu>
    class AcsrGlobalPlanner{
    public:
        static constexpr int nx_ = nx;
        static constexpr int nu_ = nu;


        AcsrGlobalPlanner()=default;


        virtual std::pair<bool,std::tuple<DM,DM,DM>> make_plan(const param_t & init, double length, int N , bool print) = 0;

        virtual std::string save(DbManager &db_manager, DM& dt_array,DM& x,DM& u,const std::vector<std::string>& x_headers={},const std::vector<std::string>& u_headers={}){
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
                os<<double(dt_array(0,i));

                for(auto j=0;j<nx_;++j){
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


    template<class TireModel = PacejkaSimpleModel>
    class BicycleDynamicsTwoBrakeOptimizer : public AcsrGlobalPlanner<7,4>{
    public:
        BicycleDynamicsTwoBrakeOptimizer() = default;

        explicit BicycleDynamicsTwoBrakeOptimizer(std::shared_ptr<Path> path_ptr,
                                         double path_width, const json & params)/*,
                                         int optimization_resolution = 100,
                                         std::shared_ptr<valid_checker_t> valid_checker=nullptr)*/
                :path_width_(path_width),path_ptr_(path_ptr){

            option_["max_iter"] = 30000;
            option_["tol"]=1e-6;
            option_["linear_solver"]="ma57";

            front_tire_model_ = std::make_shared<TireModel>(params.at("tire").at("front"));
            rear_tire_model_ = std::make_shared<TireModel>(params.at("tire").at("rear"));


            cd_ =params.at("throttle").at("Cd");
            cm0_ = params.at("throttle").at("Cm0");
            cm1_ = params.at("throttle").at("Cm1");
            cm2_ = params.at("throttle").at("Cm2");
            cbf_ = params.at("throttle").at("Cbf");
            cbr_ = params.at("throttle").at("Cbr");

            lf_ = params.at("lf");
            lr_ = params.at("lr");
            mass_ = params.at("mass");
            Iz_ = params.at("Iz");

            //const auto track_width = path_->get_width();
            v_min_ = params.at("constraint").at("v_min");
            v_max_ = params.at("constraint").at("v_max");
            d_min_ = params.at("constraint").at("d_min");
            d_max_ = params.at("constraint").at("d_max");
            delta_min_ = params.at("constraint").at("delta_min");
            delta_max_ = params.at("constraint").at("delta_max");
            delta_dot_min_ = params.at("constraint").at("delta_dot_min");
            delta_dot_max_ = params.at("constraint").at("delta_dot_max");
        }

        std::pair<bool,std::tuple<DM,DM,DM>> make_plan(const param_t & init, double length =-1, int N =1000, bool print = true){
        //std::pair<bool,std::tuple<DM,DM,DM>> make_plan(double n0 =0, double v0=0.15,bool print = true){
            auto all = Slice();
            auto _0_N = Slice(0,N);
            auto _N_1 = Slice(1,N+1);

            Dict casadi_option;
            casadi_option["print_time"]=print;
            if(print){
                option_["print_level"]=5;
            }else{
                option_["print_level"]=1;
            }

            const auto Fz = casadi::DM::ones(2) * mass_ / 2 * 9.81;

            double tau0,s0,st,phi0,v0,n0;
            if(init.find("t")!=init.end()){
                tau0 = init.at("t");
                s0 = double(path_ptr_->t_to_s_lookup(DM{tau0})[0]);
            }else{
                s0 = init.at("s");
                tau0 = double(path_ptr_->s_to_t_lookup(DM{s0})[0]);
            }

            if(init.find("phi")!=init.end()){
                phi0=init.at("phi");
            }else{
                phi0 = double(path_ptr_->f_phi(DM{tau0})[0]);
            }

            if(init.find("v")!=init.end()){
                v0=init.at("v");
            }else{
                v0 = 0.15;
            }

            if(init.find("n")!=init.end()){
                n0=init.at("n");
            }else{
                n0 = 0.0;
            }
            auto X0 = casadi::DM::vertcat({tau0, n0, phi0, v0, 0, 0, 0});
            std::cout<<"X0: "<<X0<<std::endl;

            if(length<=0){
                st = path_ptr_->get_max_length();
            }else{
                st = s0+length;
                st=std::min(st,path_ptr_->get_max_length());
            }

            const auto s_array = casadi::DM::linspace(s0,st,N+1).T();
            const auto tau_array = path_ptr_->s_to_t_lookup(s_array)[0];

            const auto tau_array_mid = path_ptr_->s_to_t_lookup((s_array(0,_0_N)+s_array(0,_N_1))*0.5);
            const auto kappa_array_mid = path_ptr_->f_kappa(tau_array_mid)[0];
            const auto tangent_vec_array_mid= path_ptr_->f_tangent_vec(tau_array_mid)[0];

            const auto phi_array_mid = DM::atan2(tangent_vec_array_mid(1,all),tangent_vec_array_mid(0,all));
            const auto tangent_vec_norm_mid = DM::sqrt((tangent_vec_array_mid(0,all)*tangent_vec_array_mid(0,all)+tangent_vec_array_mid(1,all)*tangent_vec_array_mid(1,all)));



            casadi::Opti opti;
            auto X = opti.variable(nx_,N+1);
            auto U = opti.variable(nu_, N);
            auto dt_sym_array = opti.variable(1,N);

            auto n_sym_array = (X(IDX_X_n,_0_N)+X(IDX_X_n,_N_1))*0.5;
            auto phi_sym_array = (X(IDX_X_phi, _0_N)+X(IDX_X_phi, _N_1))*0.5;
            auto vx_sym_array = (X(IDX_X_vx,_0_N)+X(IDX_X_vx,_N_1))*0.5;
            auto vy_sym_array = (X(IDX_X_vy,_0_N)+X(IDX_X_vy,_N_1))*0.5;
            auto omega_sym_array = (X(IDX_X_omega,_0_N)+X(IDX_X_omega,_N_1))*0.5;
            auto delta_sym_array = (X(IDX_X_delta,_0_N)+X(IDX_X_delta,_N_1))*0.5;

            auto delta_dot_sym_array = U(IDX_U_Ddelta,all);
            auto throttle_sym_array = U(IDX_U_Throttle,all);
            auto front_brake_sym_array = U(IDX_U_Fb,all);
            auto rear_brake_sym_array = U(IDX_U_Rb,all);

            //auto n_sym_array = X(1,Slice());
            //auto n_obj = (casadi::MX::atan(5 * ( n_sym_array*n_sym_array - track->get_width()*track->get_width() / 4) ) + casadi::pi / 2) * 12;

            auto dphi_c_sym_array = phi_sym_array - phi_array_mid;
            auto fx_f_sym_array = cd_ * throttle_sym_array - cm0_ - cm1_ * vx_sym_array * vx_sym_array - cm2_ * vy_sym_array * vy_sym_array - cbf_ * front_brake_sym_array;
            auto fx_r_sym_array = cd_ * throttle_sym_array - cm0_ - cm1_ * vx_sym_array * vx_sym_array - cm2_ * vy_sym_array * vy_sym_array - cbr_ * rear_brake_sym_array;

            auto alpha_f = MX::atan2(omega_sym_array * lf_ + vy_sym_array, -vx_sym_array) + delta_sym_array;
            auto alpha_r = MX::atan2(omega_sym_array * lr_ - vy_sym_array, vx_sym_array);

            auto fy_f_sym_array = front_tire_model_->get_lateral_force(alpha_f, double(Fz(0)));
            auto fy_r_sym_array = rear_tire_model_->get_lateral_force(alpha_r, double(Fz(1)));

            //auto n_low_obj = MX::atan(10*(MX(road_constraint_lower_bound_)-X(IDX_X_n,all)));
            //auto n_upper_obj = MX::atan(10*(X(IDX_X_n,all)-MX(road_constraint_upper_bound_)));
            auto n_obj = MX::exp(10*(X(IDX_X_n,Slice())/(path_width_/2)-1)) + MX::exp(10*(X(IDX_X_n,Slice())/(-path_width_/2)-1));

            opti.minimize(100*MX::sum2(dt_sym_array) + 10*MX::dot(delta_dot_sym_array,delta_dot_sym_array) + MX::sum2(n_obj));
            //dynamics
            //opti.subject_to(X(all, Slice(1,N+1)) == X(all, _N_1) + X_dot);
            opti.subject_to(X(0, _N_1) == X(0, _0_N) + dt_sym_array * (vx_sym_array * MX::cos(dphi_c_sym_array) - vy_sym_array * MX::sin(dphi_c_sym_array))/(tangent_vec_norm_mid*(1-n_sym_array*kappa_array_mid)));
            opti.subject_to(X(1, _N_1) == X(1, _0_N) + dt_sym_array * (vx_sym_array * MX::sin(dphi_c_sym_array) + vy_sym_array* MX::cos(dphi_c_sym_array)));
            opti.subject_to(X(2, _N_1) == X(2, _0_N) + dt_sym_array * omega_sym_array);
            opti.subject_to(X(3, _N_1) == X(3, _0_N) + dt_sym_array * (fx_r_sym_array + fx_f_sym_array * MX::cos(delta_sym_array) - fy_f_sym_array * MX::sin(delta_sym_array) + mass_ * vy_sym_array* omega_sym_array)/ mass_);
            opti.subject_to(X(4, _N_1) == X(4, _0_N) + dt_sym_array * (fy_r_sym_array + fx_f_sym_array * MX::sin(delta_sym_array) + fy_f_sym_array * MX::cos(delta_sym_array) - mass_ * vx_sym_array * omega_sym_array)/ mass_);
            opti.subject_to(X(5, _N_1) == X(5, _0_N) + dt_sym_array * (fy_f_sym_array * lf_ * MX::cos(delta_sym_array) + fx_f_sym_array * lf_ * MX::sin(delta_sym_array) - fy_r_sym_array * lr_)/ Iz_);
            opti.subject_to(X(6, _N_1) == X(6, _0_N) + dt_sym_array * delta_dot_sym_array);
            //inital conditions

            opti.subject_to(X(0, all) == tau_array);
            opti.subject_to(X(all, 0) == X0);
            opti.subject_to(dt_sym_array >0);
            opti.subject_to(X(IDX_X_vx, all) >0);

            //state boundary
            opti.subject_to(opti.bounded(delta_min_, X(IDX_X_delta,all), delta_max_));

            //control boundary
            opti.subject_to(opti.bounded(delta_dot_min_, delta_dot_sym_array, delta_dot_max_));
            opti.subject_to(opti.bounded(0, throttle_sym_array, d_max_));
            opti.subject_to(opti.bounded(0, front_brake_sym_array, 1));
            opti.subject_to(opti.bounded(0, rear_brake_sym_array, 1));

            auto X_guess = casadi::DM::zeros(nx_,N+1);
            X_guess(Slice(), 0) = X0;
            //X_guess(IDX_X_t,all) = tau_array;

            X_guess(IDX_X_phi,0) = phi0;
            const auto phi_array = path_ptr_->f_phi(tau_array)[0];
            for(auto i=1;i<N+1;++i){
                X_guess(IDX_X_phi,i) = phi_array(i);
                if(double(phi_array(i))-double(X_guess(IDX_X_phi,i-1))<-pi){
                    X_guess(IDX_X_phi,i)=X_guess(IDX_X_phi,i)+2*pi;
                }else if(double(phi_array(i))-double(X_guess(IDX_X_phi,i-1))>pi){
                    X_guess(IDX_X_phi,i)=X_guess(IDX_X_phi,i)-2*pi;
                }
            }

            X_guess(IDX_X_vx,all)=v0;
            X_guess(IDX_X_t,all) = tau_array;

            opti.set_initial(X, X_guess);
            opti.set_initial(dt_sym_array,(s_array(_N_1)-s_array(_0_N))/v0);
            opti.solver("ipopt", casadi_option, option_);
            try{
                auto sol = opti.solve();
                auto dt_array = sol.value(dt_sym_array);
                auto sol_x = sol.value(X);
                auto sol_u = sol.value(U);
                return std::make_pair(true,std::make_tuple(dt_array,sol_x,sol_u));
            }
            catch (CasadiException& e){
                std::cout<<e.what()<<std::endl;
                std::cout<<"Solve Optimal Problem Fails\n";
                return std::make_pair(false,std::make_tuple(DM{},DM{},DM{}));;
            }
        }

    public:
        constexpr static unsigned IDX_X_t = 0;
        constexpr static unsigned IDX_X_n = 1;
        constexpr static unsigned IDX_X_phi = 2;
        constexpr static unsigned IDX_X_vx = 3;
        constexpr static unsigned IDX_X_vy= 4;
        constexpr static unsigned IDX_X_omega = 5;
        constexpr static unsigned IDX_X_delta = 6;

        constexpr static unsigned IDX_U_Ddelta = 0;
        constexpr static unsigned IDX_U_Throttle = 1;
        constexpr static unsigned IDX_U_Fb = 2;
        constexpr static unsigned IDX_U_Rb = 3;

    private:
        std::shared_ptr<Path> path_ptr_;
        double path_width_{};
        std::shared_ptr<TireModel> front_tire_model_,rear_tire_model_;
        //param_t vehicle_params_,front_tire_params_,rear_tire_params_,track_params_;
        //std::string database_file_ ="../output/global_planner.db";
        //std::string datatable_name_ = std::string();
        //std::string current_datatable_name_ = std::string();
        Dict option_;
        //bool save_to_database_=false;
        //DM road_constraint_upper_bound_,road_constraint_lower_bound_;
        //const int N_{};

        //std::shared_ptr<valid_checker_t> valid_checker_;

        //constexpr static int nx = 7;
        //constexpr static int nu = 4;
        double cd_,cm0_,cm1_,cm2_,cbf_,cbr_;
        double lf_,lr_,mass_,Iz_;
        double delta_min_,delta_max_,delta_dot_min_,delta_dot_max_,v_max_,v_min_,d_max_,d_min_;

    };





    class BicycleKinematicOptimizer:public AcsrGlobalPlanner<4,2>{
    public:
        BicycleKinematicOptimizer() = default;

        BicycleKinematicOptimizer(std::shared_ptr<Path> path_ptr,
                                           double path_width,
                                           const json & params)
                :path_width_(path_width),path_ptr_(path_ptr){

            option_["max_iter"] = 3000;
            option_["tol"]=1e-9;
            option_["linear_solver"]="ma57";


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
        }

        std::pair<bool,std::tuple<DM,DM,DM>> make_plan(const param_t & init, double length =-1, int N =100, bool print = true){


            auto all = Slice();
            auto _0_N = Slice(0,N);
            auto _N_1 = Slice(1,N+1);
            Dict casadi_option;
            casadi_option["print_time"]=print;
            if(print){
                option_["print_level"]=5;
            }else{
                option_["print_level"]=1;
            }

            double tau0,s0,st,phi0,v0;
            if(init.find("t")!=init.end()){
                tau0 = init.at("t");
                s0 = double(path_ptr_->t_to_s_lookup(DM{tau0})[0]);
            }else{
                s0 = init.at("s");
                tau0 = double(path_ptr_->s_to_t_lookup(DM{s0})[0]);
            }

            if(init.find("phi")!=init.end()){
                phi0=init.at("phi");
            }else{
                phi0 = double(path_ptr_->f_phi(DM{tau0})[0]);
            }

            if(init.find("v")!=init.end()){
                v0=init.at("v");
            }else{
                v0 = 0.15;
            }

            if(length<=0){
                st = path_ptr_->get_max_length();
            }else{
                st = s0+length;
                st=std::min(st,path_ptr_->get_max_length());
            }

            const auto s_array = casadi::DM::linspace(s0,st,N+1).T();
            const auto tau_array = path_ptr_->s_to_t_lookup(s_array)[0];

            //auto tau_array_mid = path_ptr_->s_to_t_lookup((s_val(0,_0_N)+s_val(0,_N_1))*0.5)[0];

            auto kappa_array = path_ptr_->f_kappa(tau_array)[0];
            auto tangent_vec_array= path_ptr_->f_tangent_vec(tau_array)[0];
            auto phi_array = path_ptr_->f_phi(tau_array)[0];
            auto sin_phi_array = DM::sin(phi_array(0,_0_N));
            auto cos_phi_array = DM::cos(phi_array(0,_0_N));
            //auto phi_array = DM::atan2(tangent_vec_array(1,all),tangent_vec_array(0,all));
            //auto phi_array_mid = phi_array_raw+4*pi;

            auto tangent_vec_norm = DM::sqrt((tangent_vec_array(0,all)*tangent_vec_array(0,all)+tangent_vec_array(1,all)*tangent_vec_array(1,all)));

            casadi::Opti opti;
            auto X = opti.variable(nx_,N+1);
            auto U = opti.variable(nu_, N);
            auto dt_sym_array = opti.variable(1,N);

            auto n_sym_array = X(IDX_X_n,_0_N);
            auto phi_sym_array = X(IDX_X_phi, _0_N);


            auto vx_sym_array = (X(IDX_X_vx,_0_N)+X(IDX_X_vx,_N_1))*0.5;

            auto delta_sym_array = U(IDX_U_delta,all);
            auto d_sym_array = U(IDX_U_a,all);



            auto X0 = casadi::DM::vertcat({tau0, 0, phi0, v0});
            auto dphi_c_sym_array =  phi_sym_array - phi_array(0,_0_N);

            auto n_obj = MX::exp(10*(X(IDX_X_n,Slice())/(path_width_/2)-1)) + MX::exp(10*(X(IDX_X_n,Slice())/(-path_width_/2)-1));
            //opti.minimize(MX::sum2(dt_sym_array));
            opti.minimize(1000*MX::sum2(dt_sym_array)+ MX::sum2(n_obj));

            //dynamics
            opti.subject_to(X(IDX_X_t, _N_1) == X(IDX_X_t, _0_N) + dt_sym_array * (vx_sym_array * MX::cos(dphi_c_sym_array))/(tangent_vec_norm(0,_0_N)*(1-n_sym_array*kappa_array(0,_0_N))));
            opti.subject_to(X(IDX_X_n, _N_1) == X(IDX_X_n, _0_N) + dt_sym_array * (vx_sym_array * MX::sin(dphi_c_sym_array)));
            opti.subject_to(X(IDX_X_phi, _N_1) == X(IDX_X_phi, _0_N) + dt_sym_array * vx_sym_array * MX::tan(delta_sym_array)/wheel_base_);
            //opti.subject_to(X(IDX_X_t, _N_1) == X(IDX_X_t, _0_N) + dt_sym_array * (vx_sym_array * (MX::cos(phi_sym_array)*cos_phi_array+MX::sin(phi_sym_array)*sin_phi_array))/(tangent_vec_norm(0,_0_N)*(1-n_sym_array*kappa_array(0,_0_N))));
            //opti.subject_to(X(IDX_X_n, _N_1) == X(IDX_X_n, _0_N) + dt_sym_array * (vx_sym_array * (MX::sin(phi_sym_array)*cos_phi_array-MX::cos(phi_sym_array)*sin_phi_array)));

            //opti.subject_to(X(IDX_X_phi, _N_1) == X(IDX_X_phi, _0_N) + dt_sym_array * vx_sym_array * delta_sym_array/wheel_base_);
            opti.subject_to(X(IDX_X_vx, _N_1) == X(IDX_X_vx, _0_N) + dt_sym_array * (7.9*d_sym_array-0.3*vx_sym_array));

            //inital conditions
            opti.subject_to(X(0, all) == tau_array);
            opti.subject_to(X(all, 0) == X0);
            opti.subject_to(dt_sym_array >0.0001);

            //state boundary
            //opti.subject_to(opti.bounded(-path_width_/2, X(IDX_X_n,_0_N), path_width_/2));
            opti.subject_to(vx_sym_array>=v_min_);


            //control boundary
            opti.subject_to(opti.bounded(d_min_, d_sym_array, d_max_));
            opti.subject_to(opti.bounded(delta_min_, delta_sym_array, delta_max_));


            auto X_guess = casadi::DM::zeros(nx_,N+1);
            auto U_guess = casadi::DM::zeros(nu_,N);



            X_guess(Slice(), 0) = X0;
            //X_guess(IDX_X_t,all) = tau_array;

            X_guess(IDX_X_phi,0) = phi0;
            for(auto i=1;i<N+1;++i){
                X_guess(IDX_X_phi,i) = phi_array(i);
                if(double(phi_array(i))-double(X_guess(IDX_X_phi,i-1))<-pi){
                    X_guess(IDX_X_phi,i)=X_guess(IDX_X_phi,i)+2*pi;
                }else if(double(phi_array(i))-double(X_guess(IDX_X_phi,i-1))>pi){
                    X_guess(IDX_X_phi,i)=X_guess(IDX_X_phi,i)-2*pi;
                }
            }


            X_guess(IDX_X_vx,all)=v0;

            opti.set_initial(X, X_guess);
            opti.set_initial(U, U_guess);
            opti.set_initial(dt_sym_array,(s_array(_N_1)-s_array(_0_N))/v0);

            opti.solver("ipopt", casadi_option, option_);
            try{
                auto sol = opti.solve();
                auto dt_array = sol.value(dt_sym_array);
                auto sol_x = sol.value(X);
                auto sol_u = sol.value(U);
                return std::make_pair(true,std::make_tuple(dt_array,sol_x,sol_u));
            }
            catch (CasadiException& e){
                std::cout<<e.what()<<std::endl;
                std::cout<<"Solve Optimal Problem Fails\n";
                return std::make_pair(false,std::make_tuple(DM{},DM{},DM{}));
            }

        }

    public:
        constexpr static unsigned IDX_X_t = 0;
        constexpr static unsigned IDX_X_n = 1;
        constexpr static unsigned IDX_X_phi = 2;
        constexpr static unsigned IDX_X_vx = 3;
        constexpr static unsigned IDX_U_delta = 0;
        constexpr static unsigned IDX_U_a = 1;

    private:
        std::shared_ptr<Path> path_ptr_;
        double path_width_{};
        Dict option_;
        double delta_min_,delta_max_,v_min_,v_max_,d_min_,d_max_,wheel_base_;

    };



    /*
    class BicycleDynamicsOneBrakeOptimizer{
    public:
        BicycleDynamicsOneBrakeOptimizer() = default;
        BicycleDynamicsOneBrakeOptimizer(const DM& waypoints,double max_width,const std::string& vehicle_file,const std::string& front_tire_file,const std::string& rear_tire_file){
            if(!std::filesystem::exists(vehicle_file)){
                std::cout<<vehicle_file<<" does not exist\n";
                return;
            }
            if(!std::filesystem::exists(front_tire_file)){
                std::cout<<front_tire_file<<" does not exist\n";
                return;
            }
            if(!std::filesystem::exists(rear_tire_file)){
                std::cout<<front_tire_file<<" does not exist\n";
                return;
            }

            //load params
            auto vehicle_yaml = YAML::LoadFile(vehicle_file);
            std::cout<<"load vehicle config file.. Total node: "<<vehicle_yaml.size()<<std::endl;
            for(auto it = vehicle_yaml.begin();it!=vehicle_yaml.end();++it){
                vehicle_params[it->first.as<std::string>()]=it->second.as<double>();
            }

            auto front_tire_yaml = YAML::LoadFile(front_tire_file);
            std::cout<<"load front tire config file.. Total node: "<<front_tire_yaml.size()<<std::endl;
            for(auto it = front_tire_yaml.begin();it!=front_tire_yaml.end();++it){
                front_tire_params[it->first.as<std::string>()]=it->second.as<double>();
            }

            if(front_tire_file.compare(rear_tire_file)==0){
                rear_tire_params = front_tire_params;
            }else {
                auto rear_tire_yaml = YAML::LoadFile(rear_tire_file);
                std::cout<<"load rear tire config file.. Total node: "<<rear_tire_yaml.size()<<std::endl;
                for (auto it = rear_tire_yaml.begin(); it != rear_tire_yaml.end(); ++it) {
                    rear_tire_params[it->first.as<std::string>()] = it->second.as<double>();
                }
            }

            path_ = std::make_shared<Path>(waypoints);
            option["max_iter"] = 10000;
            option["tol"]=1e-6;
            option["linear_solver"]="ma57";
            //option["print_level"]=5;



        }

        void make_plan(double start_s = 0, double end_s = -1,double n0 =0, double v0=0.15,int N = 100,bool print = true, bool save_to_database = true){
            auto all = Slice();
            auto _0_N = Slice(0,N);
            auto _N_1 = Slice(1,N+1);
            Dict casadi_option;
            casadi_option["print_time"]=print;
            if(print){
                option["print_level"]=5;
            }else{
                option["print_level"]=3;
            }

            auto front_tire_model = std::make_shared<PacejkaSimpleModel>(front_tire_params);
            auto rear_tire_model = std::make_shared<PacejkaSimpleModel>(rear_tire_params);

            const auto cd =vehicle_params.at("Cd");
            const auto cm0 = vehicle_params.at("Cm0");
            const auto cm1 = vehicle_params.at("Cm1");
            const auto cm2 = vehicle_params.at("Cm2");
            const auto cbf = vehicle_params.at("Cbf");
            const auto cbr = vehicle_params.at("Cbr");
            const auto lf = vehicle_params.at("lf");
            const auto lr = vehicle_params.at("lr");
            const auto mass = vehicle_params.at("m");
            const auto Iz = vehicle_params.at("Iz");
            const auto Fz = casadi::DM::ones(2) * vehicle_params.at("m") / 2 * 9.81;

            const auto track_width = path_->get_width();
            const auto delta_min = vehicle_params.at("delta_min");
            const auto delta_max = vehicle_params.at("delta_max");
            const auto delta_dot_min = vehicle_params.at("delta_dot_min");
            const auto delta_dot_max = vehicle_params.at("delta_dot_max");

            const auto s_val = casadi::DM::linspace(start_s,end_s,N+1).T();
            const auto tau_array = path_->s_to_t_lookup(s_val)[0];
            const auto tau0 = double(tau_array(0));
            const auto phi0 = double(path_->f_phi(tau_array(0))[0]);
            const auto tau_array_mid = path_->s_to_t_lookup((s_val(0,_0_N)+s_val(0,_N_1))*0.5);

            const auto kappa_array = path_->f_kappa(tau_array_mid)[0];
            //const auto kappa_array = (kappa_array_temp(0,_N_1) + kappa_array_temp(0,Slice(1,N+1)))*0.5;
            const auto tangent_vec_array= path_->f_tangent_vec(tau_array_mid)[0];
            //const auto tangent_vec_array = (tangent_vec_array_temp(all,_N_1) + tangent_vec_array_temp(all,Slice(1,N+1)))*0.5;

            const auto phi_array = DM::atan2(tangent_vec_array(1,all),tangent_vec_array(0,all));
            const auto tangent_vec_norm = DM::sqrt((tangent_vec_array(0,all)*tangent_vec_array(0,all)+tangent_vec_array(1,all)*tangent_vec_array(1,all)));

            //tau_array.to_file("/home/acsr/Documents/data/tau_array_cpp.txt");
            //kappa_array.to_file("/home/acsr/Documents/data/kappa_array_cpp.txt");
            //tangent_vec_array.to_file("/home/acsr/Documents/data/tangent_vec_array_cpp.txt");
            //phi_array.to_file("/home/acsr/Documents/data/phi_array_cpp.txt");
            //tangent_vec_norm.to_file("/home/acsr/Documents/data/tangent_vec_norm_cpp.txt");

            casadi::Opti opti;
            auto X = opti.variable(nx,N+1);
            auto U = opti.variable(nu, N);
            //auto X_dot = opti.variable(nx, N);
            auto dt_sym_array = opti.variable(1,N);

            auto n_sym_array = (X(IDX_X_n,_0_N)+X(IDX_X_n,_N_1))*0.5;
            auto phi_sym_array = (X(IDX_X_phi, _0_N)+X(IDX_X_phi, _N_1))*0.5;
            auto vx_sym_array = (X(IDX_X_vx,_0_N)+X(IDX_X_vx,_N_1))*0.5;
            auto vy_sym_array = (X(IDX_X_vy,_0_N)+X(IDX_X_vy,_N_1))*0.5;
            auto omega_sym_array = (X(IDX_X_omega,_0_N)+X(IDX_X_omega,_N_1))*0.5;
            auto delta_sym_array = (X(IDX_X_delta,_0_N)+X(IDX_X_delta,_N_1))*0.5;

            auto delta_dot_sym_array = U(IDX_U_Ddelta,all);
            auto throttle_sym_array = U(IDX_U_Throttle,all);
            auto brake_sym_array = U(IDX_U_Brake,all);
            //auto rear_brake_sym_array = U(IDX_U_Rb,all);

            //auto n_sym_array = X(1,Slice());
            //auto n_obj = (casadi::MX::atan(5 * ( n_sym_array*n_sym_array - track->get_width()*track->get_width() / 4) ) + casadi::pi / 2) * 12;

            auto X0 = casadi::DM::vertcat({tau0, n0, phi0, v0, 0, 0, 0});
            std::cout<<"X0: "<<X0<<std::endl;

            auto dphi_c_sym_array = phi_sym_array - phi_array;
            auto fx_sym_array = cd * throttle_sym_array - cm0 - cm1 * vx_sym_array * vx_sym_array - cm2 * vy_sym_array * vy_sym_array - cbf * brake_sym_array;
            //auto fx_r_sym_array = cd * throttle_sym_array - cm0 - cm1 * vx_sym_array * vx_sym_array - cm2 * vy_sym_array * vy_sym_array - cbr * brake_sym_array;

            auto alpha_f = MX::atan2(omega_sym_array * lf + vy_sym_array, -vx_sym_array) + delta_sym_array;
            auto alpha_r = MX::atan2(omega_sym_array * lr - vy_sym_array, vx_sym_array);

            auto fy_f_sym_array = front_tire_model->get_lateral_force(alpha_f, double(Fz(0)));
            auto fy_r_sym_array = rear_tire_model->get_lateral_force(alpha_r, double(Fz(1)));

            //X_dot(0, all) = dt_sym_array*(vx_sym_array * MX::cos(dphi_c_sym_array) - vy_sym_array * MX::sin(dphi_c_sym_array))/(tangent_vec_norm*(1-n_sym_array*kappa_array));
            //X_dot(1, all) = dt_sym_array*(vx_sym_array * MX::sin(dphi_c_sym_array) + vy_sym_array* MX::cos(dphi_c_sym_array)) ;
            //X_dot(2, all) = dt_sym_array*omega_sym_array;
            //X_dot(3, all) = dt_sym_array * (fx_r_sym_array + fx_f_sym_array * MX::cos(delta_sym_array) - fy_f_sym_array * MX::sin(delta_sym_array) + mass * vy_sym_array* omega_sym_array)/ mass;
            //X_dot(4, all) = dt_sym_array * (fy_r_sym_array + fx_f_sym_array * MX::sin(delta_sym_array) + fy_f_sym_array * MX::cos(delta_sym_array) - mass * vx_sym_array * omega_sym_array)/ mass;
            //X_dot(5, all) = dt_sym_array * (fy_f_sym_array * lf * MX::cos(delta_sym_array) + fx_f_sym_array * lf * MX::sin(delta_sym_array) - fy_r_sym_array * lr)/ Iz ;
            //X_dot(6, all) = dt_sym_array * delta_dot_sym_array;



            auto n_obj = (MX::atan(5 * (n_sym_array * n_sym_array - (track_width / 2) *(track_width / 2))) + casadi::pi / 2) * 4.0;
            //opti.minimize(MX::sum2(dt_sym_array) + MX::dot(n_obj, n_obj));
            opti.minimize(MX::sum2(dt_sym_array) + 0.1*MX::dot(delta_dot_sym_array,delta_dot_sym_array) + MX::dot(n_obj,n_obj));

            //dynamics
            //opti.subject_to(X(all, Slice(1,N+1)) == X(all, _N_1) + X_dot);
            opti.subject_to(X(0, _N_1) == X(0, _0_N) + dt_sym_array * (vx_sym_array * MX::cos(dphi_c_sym_array) - vy_sym_array * MX::sin(dphi_c_sym_array))/(tangent_vec_norm*(1-n_sym_array*kappa_array)));
            opti.subject_to(X(1, _N_1) == X(1, _0_N) + dt_sym_array * (vx_sym_array * MX::sin(dphi_c_sym_array) + vy_sym_array* MX::cos(dphi_c_sym_array)));
            opti.subject_to(X(2, _N_1) == X(2, _0_N) + dt_sym_array * omega_sym_array);
            opti.subject_to(X(3, _N_1) == X(3, _0_N) + dt_sym_array * (fx_sym_array + fx_sym_array * MX::cos(delta_sym_array) - fy_f_sym_array * MX::sin(delta_sym_array) + mass * vy_sym_array* omega_sym_array)/ mass);
            opti.subject_to(X(4, _N_1) == X(4, _0_N) + dt_sym_array * (fy_r_sym_array + fx_sym_array * MX::sin(delta_sym_array) + fy_f_sym_array * MX::cos(delta_sym_array) - mass * vx_sym_array * omega_sym_array)/ mass);
            opti.subject_to(X(5, _N_1) == X(5, _0_N) + dt_sym_array * (fy_f_sym_array * lf * MX::cos(delta_sym_array) + fx_sym_array * lf * MX::sin(delta_sym_array) - fy_r_sym_array * lr)/ Iz);
            opti.subject_to(X(6, _N_1) == X(6, _0_N) + dt_sym_array * delta_dot_sym_array);
            //inital conditions

            opti.subject_to(X(0, all) == tau_array);
            opti.subject_to(X(all, 0) == X0);
            opti.subject_to(dt_sym_array >0);

            //state boundary
            opti.subject_to(opti.bounded(delta_min, X(IDX_X_delta,all), delta_max));
            //opti.subject_to(opti.bounded(-track_width/2,X(IDX_X_n,all),track_width/2));
            //opti.subject_to(opti.bounded(delta_min, delta_sym_array, delta_max));
            //opti.subject_to(opti.bounded(-track_width/2,n_sym_array,track_width/2));
            //control boundary
            //opti.subject_to(opti.bounded(delta_dot_min, delta_dot_sym_array, delta_dot_max));
            opti.subject_to(opti.bounded(0, throttle_sym_array, 1));
            opti.subject_to(opti.bounded(0, brake_sym_array, 1));

            auto X_guess = casadi::DM::zeros(nx,N+1);
            //X_guess(0,Slice()) = tau_array;
            //X_guess(2,Slice(1,N+1)) = phi_array;
            X_guess(3,all) = v0;
            //for(auto i=0;i<N+1;++i){
            //    X_guess(Slice(), i) = X0;
            //}
            opti.set_initial(X, X_guess);

            opti.solver("ipopt", casadi_option, option);
            try{
                auto sol = opti.solve();
                if(save_to_database){
                    auto dt_array = sol.value(dt_sym_array);
                    auto sol_x = sol.value(X);
                    auto sol_u = sol.value(U);
                    save(dt_array,sol_x,sol_u);
                }
            }
            catch (CasadiException& e){
                std::cout<<e.what()<<std::endl;
                std::cout<<"Solve Optimal Problem Fails\n";
                return;
            }

        }


        void set_database_file(const std::string& file_name){
            database_file =  file_name;
        }

        void set_datatable_name(const std::string& table_name){
            datatable_name = table_name;
        }

    private:
        std::shared_ptr<Path> path_;
        param_t vehicle_params,front_tire_params,rear_tire_params,track_params;
        std::string database_file ="../output/global_planner.db";
        std::string datatable_name = std::string();
        std::string current_datatable_name = std::string();
        Dict option;


        constexpr static int nx = 7;
        constexpr static int nu = 3;
        constexpr static unsigned IDX_X_t = 0;
        constexpr static unsigned IDX_X_n = 1;
        constexpr static unsigned IDX_X_phi = 2;
        constexpr static unsigned IDX_X_vx = 3;
        constexpr static unsigned IDX_X_vy= 4;
        constexpr static unsigned IDX_X_omega = 5;
        constexpr static unsigned IDX_X_delta = 6;

        constexpr static unsigned IDX_U_Ddelta = 0;
        constexpr static unsigned IDX_U_Throttle = 1;
        constexpr static unsigned IDX_U_Brake = 2;

        void save(DM& dt_array,DM& x,DM& u){


            SQLite::Database    db(database_file, SQLite::OPEN_READWRITE|SQLite::OPEN_CREATE);
            std::cout << "SQLite database file '" << db.getFilename().c_str() << "' opened successfully\n";

            current_datatable_name = datatable_name;

            if(current_datatable_name.empty()){
                const auto now = std::chrono::system_clock::now();
                time_t rawtime;
                struct tm * timeinfo;
                char buffer[40];
                time (&rawtime);
                timeinfo = localtime(&rawtime);
                strftime(buffer,sizeof(buffer),"_%d_%H_%M_%S",timeinfo);
                current_datatable_name = std::string(buffer);
            }


            {
                //delete table if exist
                std::string drop_statement_string = "DROP TABLE IF EXISTS " + current_datatable_name;
                db.exec(drop_statement_string);
                std::cout << "Save data to data table: " << current_datatable_name << "\n";

                //create table
                std::string table_statement_string = "create table if not exists "
                                                     +current_datatable_name
                                                     + "("
                                                       "id INTEGER PRIMARY KEY AUTOINCREMENT, "
                                                       "dt REAL,"
                                                       "tau REAL,"
                                                       "s REAL,"
                                                       "n REAL,"
                                                       "x REAL,"
                                                       "y REAL,"
                                                       "phi REAL,"
                                                       "vx REAL,"
                                                       "vy REAL,"
                                                       "omega REAL,"
                                                       "delta REAL,"
                                                       "delta_dot REAL,"
                                                       "throttle REAL,"
                                                       //"front_brake REAL,"
                                                       "brake REAL"
                                                       ")";
                SQLite::Statement query(db, table_statement_string);
                try {
                    query.exec();
                } catch (SQLite::Exception &e) {
                    std::cout << "Create Table "<<current_datatable_name<<" Error\n";
                    std::cout << e.what()<<std::endl;
                    return;
                }
                std::string query_string = "INSERT INTO " +current_datatable_name
                                           +" (dt,tau,s,n,x,y,phi,vx,vy,omega,delta,delta_dot,throttle,brake) "
                                            "VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?)";
                query = SQLite::Statement(db, query_string);


                auto N = x.columns()-1;
                auto dm_xy = path_->f_tn_to_xy(std::vector<DM>{x(0,Slice()),x(1,Slice())})[0];
                auto s_array = path_->t_to_s_lookup(x(0,Slice()))[0];

                for(auto i=0;i<N;++i){
                    query = SQLite::Statement(db, query_string);
                    query.bind(1,double(dt_array(0,i)));
                    query.bind(2,double(x(0,i)));
                    query.bind(3,double(s_array(0,i)));
                    query.bind(4,double(x(1,i)));
                    query.bind(5,double(dm_xy(0,i)));
                    query.bind(6,double(dm_xy(1,i)));
                    query.bind(7,double(x(2,i)));
                    query.bind(8,double(x(3,i)));
                    query.bind(9,double(x(4,i)));
                    query.bind(10,double(x(5,i)));
                    query.bind(11,double(x(6,i)));
                    query.bind(12,double(u(0,i)));
                    query.bind(13,double(u(1,i)));
                    query.bind(14,double(u(2,i)));
                    try {
                        query.exec();
                    } catch (SQLite::Exception &e) {
                        std::cout << "Insert Solution Error\n";
                        std::cout << e.what();
                        return;
                    }
                }

                query = SQLite::Statement(db, query_string);
                //query.bind(1,NULL);
                query.bind(2,double(x(0,N)));
                query.bind(3,double(s_array(0,N)));
                query.bind(4,double(x(1,N)));
                query.bind(5,double(dm_xy(0,N)));
                query.bind(6,double(dm_xy(1,N)));
                query.bind(7,double(x(2,N)));
                query.bind(8,double(x(3,N)));
                query.bind(9,double(x(4,N)));
                query.bind(10,double(x(5,N)));
                query.bind(11,double(x(6,N)));
                //query.bind(12,double(u(0,NULL)));
                //query.bind(13,double(u(1,NULL)));
                //query.bind(14,double(u(2,NULL)));
                //query.bind(15,double(u(3,NULL)));
                try {
                    query.exec();
                } catch (SQLite::Exception &e) {
                    std::cout << "Insert last dataset Error\n";
                    std::cout << e.what();
                    return;
                }

            }
            std::cout<<"save to database finished\n";


        }

    };*/
}

#endif //RACECARPLANNER_GLOBAL_PLANNER_OPTIMIZER_HPP
