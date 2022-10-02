//
// Created by acsr on 9/27/22.
//

#ifndef RACECARPLANNER_GLOBAL_PLANNER_OPTIMIZER_HPP
#define RACECARPLANNER_GLOBAL_PLANNER_OPTIMIZER_HPP

#include "path.hpp"
#include <yaml-cpp/yaml.h>
#include <filesystem>
#include "dynamics.hpp"
#include "tire_model.hpp"
#include <SQLiteCpp/SQLiteCpp.h>
#include <boost/format.hpp>
#include <filesystem>
#include <execution>

namespace acsr{
    using param_t = std::map<std::string,double>;

    template<class valid_checker_t=void>
    class BicycleDynamicsTwoBrakeOptimizer{
    public:
        BicycleDynamicsTwoBrakeOptimizer() = default;

        explicit BicycleDynamicsTwoBrakeOptimizer(std::shared_ptr<Path> path_ptr,
                                         double path_width,
                                         const param_t & vehicle_params,
                                         const param_t & front_tire_params,
                                         const param_t & rear_tire_params,
                                         int optimization_resolution = 100,
                                         std::shared_ptr<valid_checker_t> valid_checker=nullptr)
                :vehicle_params_(vehicle_params),path_width_(path_width),front_tire_params_(front_tire_params),rear_tire_params_(rear_tire_params),N_(optimization_resolution),path_(path_ptr){

            option_["max_iter"] = 30000;
            option_["tol"]=1e-6;
            option_["linear_solver"]="ma57";


            if(std::is_void<valid_checker_t>::value || valid_checker== nullptr){
                road_constraint_upper_bound_=DM{path_width/2};
                road_constraint_lower_bound_=DM{-path_width/2};
            }else{
                std::vector<size_t> index_vect(N_+1);
                std::iota(index_vect.begin(),index_vect.end(),0);

                auto s = DM::linspace(0,path_->get_max_length(),N_+1);
                s = s.T();
                auto t = path_->s_to_t_lookup(s)[0];
                auto n = DM::zeros(1,N_+1);
                auto center_line = path_->f_tn_to_xy(DMVector {t,n})[0];

                road_constraint_upper_bound_=path_width/2*DM::ones(1,N_+1);
                road_constraint_lower_bound_=-path_width/2*DM::ones(1,N_+1);

                std::for_each(std::execution::par,index_vect.begin(),index_vect.end(),[&](size_t idx){
                    std::vector<std::pair<double,DM>> result = valid_checker->template near_collides<DM>(path_width/2,double(center_line(0,idx)),double(center_line(1,idx)));
                    for(auto& r: result){
                        auto n = path_->f_xy_to_tn(DMVector{r.second,t(0,idx)})[0](0);
                        if(double(n(0))<0 && -r.first>double(road_constraint_lower_bound_(0,idx))){
                            road_constraint_lower_bound_(0,idx)=-r.first;
                        }else if(double(n(0))>0 && r.first<double(road_constraint_upper_bound_(0,idx))){
                            road_constraint_upper_bound_(0,idx) = r.first;
                        }
                    }
                });

                if(path_width/2 - double(DM::mmin(road_constraint_upper_bound_))<0.1 && path_width/2 + double(DM::mmax(road_constraint_upper_bound_))<0.1){
                    road_constraint_upper_bound_=DM{path_width/2};
                    road_constraint_lower_bound_=DM{-path_width/2};
                }

            }

        }

        std::pair<bool,std::tuple<DM,DM,DM>> make_plan(double n0 =0, double v0=0.15,bool print = true){
            auto all = Slice();
            auto _0_N = Slice(0,N_);
            auto _N_1 = Slice(1,N_+1);
            Dict casadi_option;
            casadi_option["print_time"]=print;
            if(print){
                option_["print_level"]=5;
            }else{
                option_["print_level"]=1;
            }

            auto front_tire_model = std::make_shared<PacejkaSimpleModel>(front_tire_params_);
            auto rear_tire_model = std::make_shared<PacejkaSimpleModel>(rear_tire_params_);

            const auto cd =vehicle_params_.at("Cd");
            const auto cm0 = vehicle_params_.at("Cm0");
            const auto cm1 = vehicle_params_.at("Cm1");
            const auto cm2 = vehicle_params_.at("Cm2");
            const auto cbf = vehicle_params_.at("Cbf");
            const auto cbr = vehicle_params_.at("Cbr");
            const auto lf = vehicle_params_.at("lf");
            const auto lr = vehicle_params_.at("lr");
            const auto mass = vehicle_params_.at("m");
            const auto Iz = vehicle_params_.at("Iz");
            const auto Fz = casadi::DM::ones(2) * vehicle_params_.at("m") / 2 * 9.81;

            //const auto track_width = path_->get_width();
            const auto delta_min = vehicle_params_.at("delta_min");
            const auto delta_max = vehicle_params_.at("delta_max");
            const auto delta_dot_min = vehicle_params_.at("delta_dot_min");
            const auto delta_dot_max = vehicle_params_.at("delta_dot_max");

            const auto s_val = casadi::DM::linspace(0,path_->get_max_length(),N_+1).T();
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
            auto X = opti.variable(nx,N_+1);
            auto U = opti.variable(nu, N_);
            //auto X_dot = opti.variable(nx, N);
            auto dt_sym_array = opti.variable(1,N_);

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

            auto X0 = casadi::DM::vertcat({tau0, n0, phi0, v0, 0, 0, 0});
            std::cout<<"X0: "<<X0<<std::endl;

            auto dphi_c_sym_array = phi_sym_array - phi_array;
            auto fx_f_sym_array = cd * throttle_sym_array - cm0 - cm1 * vx_sym_array * vx_sym_array - cm2 * vy_sym_array * vy_sym_array - cbf * front_brake_sym_array;
            auto fx_r_sym_array = cd * throttle_sym_array - cm0 - cm1 * vx_sym_array * vx_sym_array - cm2 * vy_sym_array * vy_sym_array - cbr * rear_brake_sym_array;

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



            //auto n_low_obj = MX::exp(MX(road_constraint_lower_bound_)-X(IDX_X_n,all));
            //auto n_upper_obj = MX::exp(X(IDX_X_n,all)-MX(road_constraint_upper_bound_));
            auto n_low_obj = MX::atan(10*(MX(road_constraint_lower_bound_)-X(IDX_X_n,all)));
            auto n_upper_obj = MX::atan(10*(X(IDX_X_n,all)-MX(road_constraint_upper_bound_)));

            //auto n_obj = MX::atan(5 * (n_sym_array * n_sym_array - (path_width_ / 2) *(path_width_ / 2))) + casadi::pi / 2;
            //opti.minimize(MX::sum2(dt_sym_array) + MX::dot(n_obj, n_obj));
            //opti.minimize(MX::sum2(dt_sym_array) + MX::dot(delta_dot_sym_array,delta_dot_sym_array) + 15.0*MX::dot(n_obj,n_obj));
            //opti.minimize(MX::sum2(dt_sym_array) + MX::dot(delta_dot_sym_array,delta_dot_sym_array));
            opti.minimize(10*MX::sum2(dt_sym_array) + 10*MX::dot(delta_dot_sym_array,delta_dot_sym_array) + MX::sum2(n_low_obj) + MX::sum2(n_upper_obj));
            //dynamics
            //opti.subject_to(X(all, Slice(1,N+1)) == X(all, _N_1) + X_dot);
            opti.subject_to(X(0, _N_1) == X(0, _0_N) + dt_sym_array * (vx_sym_array * MX::cos(dphi_c_sym_array) - vy_sym_array * MX::sin(dphi_c_sym_array))/(tangent_vec_norm*(1-n_sym_array*kappa_array)));
            opti.subject_to(X(1, _N_1) == X(1, _0_N) + dt_sym_array * (vx_sym_array * MX::sin(dphi_c_sym_array) + vy_sym_array* MX::cos(dphi_c_sym_array)));
            opti.subject_to(X(2, _N_1) == X(2, _0_N) + dt_sym_array * omega_sym_array);
            opti.subject_to(X(3, _N_1) == X(3, _0_N) + dt_sym_array * (fx_r_sym_array + fx_f_sym_array * MX::cos(delta_sym_array) - fy_f_sym_array * MX::sin(delta_sym_array) + mass * vy_sym_array* omega_sym_array)/ mass);
            opti.subject_to(X(4, _N_1) == X(4, _0_N) + dt_sym_array * (fy_r_sym_array + fx_f_sym_array * MX::sin(delta_sym_array) + fy_f_sym_array * MX::cos(delta_sym_array) - mass * vx_sym_array * omega_sym_array)/ mass);
            opti.subject_to(X(5, _N_1) == X(5, _0_N) + dt_sym_array * (fy_f_sym_array * lf * MX::cos(delta_sym_array) + fx_f_sym_array * lf * MX::sin(delta_sym_array) - fy_r_sym_array * lr)/ Iz);
            opti.subject_to(X(6, _N_1) == X(6, _0_N) + dt_sym_array * delta_dot_sym_array);
            //inital conditions

            opti.subject_to(X(0, all) == tau_array);
            opti.subject_to(X(all, 0) == X0);
            opti.subject_to(dt_sym_array >0);
            opti.subject_to(X(IDX_X_vx, all) >0);
            //state boundary
            opti.subject_to(opti.bounded(delta_min, X(IDX_X_delta,all), delta_max));
            //opti.subject_to(opti.bounded(road_constraint_lower_bound_,X(IDX_X_n,all),road_constraint_upper_bound_));
            //opti.subject_to(opti.bounded(delta_min, delta_sym_array, delta_max));
            //opti.subject_to(opti.bounded(-track_width/2,n_sym_array,track_width/2));
            //control boundary
            opti.subject_to(opti.bounded(delta_dot_min, delta_dot_sym_array, delta_dot_max));
            opti.subject_to(opti.bounded(0, throttle_sym_array, 1));
            opti.subject_to(opti.bounded(0, front_brake_sym_array, 1));
            opti.subject_to(opti.bounded(0, rear_brake_sym_array, 1));

            auto X_guess = casadi::DM::zeros(nx,N_+1);
            //X_guess(0,Slice()) = tau_array;
            //X_guess(2,Slice(1,N+1)) = phi_array;
            //X_guess(3,all) = v0;
            for(auto i=0;i<N_+1;++i){
                X_guess(Slice(), i) = X0;
            }
            opti.set_initial(X, X_guess);

            opti.solver("ipopt", casadi_option, option_);
            try{
                auto sol = opti.solve();
                auto dt_array = sol.value(dt_sym_array);
                auto sol_x = sol.value(X);
                auto sol_u = sol.value(U);
                if(save_to_database_){
                    save(dt_array,sol_x,sol_u);
                }
                return std::make_pair(true,std::make_tuple(dt_array,sol_x,sol_u));
            }
            catch (CasadiException& e){
                std::cout<<e.what()<<std::endl;
                std::cout<<"Solve Optimal Problem Fails\n";
                return std::make_pair(false,std::make_tuple(DM{},DM{},DM{}));;
            }

        }

        std::pair<DM,DM> boundary_xy(){
            auto s = DM::linspace(0,path_->get_max_length(),N_+1);
            s = s.T();
            auto t = path_->s_to_t_lookup(s)[0];
            //auto n = DM::zeros(1,N_+1);
            auto outer_line = path_->f_tn_to_xy(DMVector {t,road_constraint_lower_bound_})[0];
            auto inner_line = path_->f_tn_to_xy(DMVector {t,road_constraint_upper_bound_})[0];
            return std::make_pair(outer_line,inner_line);
        }

        std::pair<DM,DM> boundary_n(){
            return std::make_pair(road_constraint_lower_bound_,road_constraint_upper_bound_);
        }

        void set_database_file(const std::string& file_name){
            database_file_ =  file_name;
        }

        void set_datatable_name(const std::string& table_name){
            datatable_name_ = table_name;
        }


        /*
        void plot_trajectory(const std::shared_ptr<Track>& plot_track = nullptr,const std::string& table_name=std::string()){
            namespace plt = matplotlibcpp;
            std::string table_;
            std::shared_ptr<Track> track_;
            if(current_datatable_name.empty() && table_name.empty()){
                std::cout<<"no data find\n";
                return;
            }
            if(!table_name.empty()){
                table_ = table_name;
            }else{
                table_=current_datatable_name;
            }

            if(plot_track== nullptr && track==nullptr){
                std::cout<<"no track data\n";
                return;
            }

            if(plot_track)track_=plot_track;
            else track_ = track;

            SQLite::Database    db(database_file, SQLite::OPEN_READONLY);
            std::cout << "SQLite database file '" << db.getFilename().c_str() << "' opened successfully\n";
            //std::string count_str = db.execAndGet("SELECT COUNT(*) FROM " + table_);
            //int count = std::stoi(count_str);

            std::vector<double> tau_vec;
            std::vector<double> s_vec;
            std::vector<double> x_vec;
            std::vector<double> y_vec;
            std::vector<double> phi_vec;
            std::vector<double> vx_vec;
            std::vector<double> vy_vec;
            std::vector<double> n_vec;
            std::vector<double> dt_vec;
            std::vector<double> delta_vec;
            std::vector<double> throttle_vec;
            std::vector<double> br_vec;
            std::vector<double> bf_vec;

            SQLite::Statement query(db, "SELECT dt,tau,s,n,x,y,phi,vx,vy,delta,throttle,front_brake,rear_brake FROM "+table_);
            while (query.executeStep())
            {
                dt_vec.push_back(query.getColumn(0));
                tau_vec.push_back(query.getColumn(1));
                s_vec.push_back(query.getColumn(2));
                n_vec.push_back(query.getColumn(3));
                x_vec.push_back(query.getColumn(4));
                y_vec.push_back(query.getColumn(5));
                phi_vec.push_back(query.getColumn(6));
                vx_vec.push_back(query.getColumn(7));
                vy_vec.push_back(query.getColumn(8));
                delta_vec.push_back(query.getColumn(9));
                throttle_vec.push_back(query.getColumn(10));
                bf_vec.push_back(query.getColumn(11));
                br_vec.push_back(query.getColumn(12));
            }

            auto tau_vec_for_track = tau_vec;
            int size = tau_vec.size()/10;
            double dtau = (tau_vec.back()-tau_vec.front())/(tau_vec.size()-1);
            for(int i=0;i<size;++i){
                tau_vec_for_track.insert(tau_vec_for_track.begin(),tau_vec_for_track.front()-dtau);
                tau_vec_for_track.push_back(tau_vec_for_track.back()+dtau);
            }
            std::transform(tau_vec_for_track.begin(),tau_vec_for_track.end(),tau_vec_for_track.begin(),[&track_](double v){
                if(v<0)v+=track_->get_max_tau();
                return v;
            });

            DM center_line_tau = (DM(tau_vec_for_track)).T();
            DM zero_vec = DM::zeros(1,tau_vec_for_track.size());
            DM one_vec = DM::ones(1,tau_vec_for_track.size());
            auto center_line = track_->f_tn_to_xy(DMVector {center_line_tau,zero_vec})[0];
            auto inner_line = track_->f_tn_to_xy(DMVector {center_line_tau,track_->get_width()/2*one_vec})[0];
            auto outer_line = track_->f_tn_to_xy(DMVector {center_line_tau,-track_->get_width()/2*one_vec})[0];

            auto center_x_dm = center_line(0,Slice());
            auto center_y_dm = center_line(1,Slice());
            auto inner_x_dm = inner_line(0,Slice());
            auto inner_y_dm = inner_line(1,Slice());
            auto outer_x_dm = outer_line(0,Slice());
            auto outer_y_dm = outer_line(1,Slice());
            plt::figure(1);
            plt::plot(std::vector<double>{center_x_dm->begin(),center_x_dm->end()},
                      std::vector<double>{center_y_dm->begin(),center_y_dm->end()},
                      "y--",
                      std::vector<double>{inner_x_dm->begin(),inner_x_dm->end()},
                      std::vector<double>{inner_y_dm->begin(),inner_y_dm->end()},
                      "k-",
                      std::vector<double>{outer_x_dm->begin(),outer_x_dm->end()},
                      std::vector<double>{outer_y_dm->begin(),outer_y_dm->end()},
                      "k-",
                      x_vec,
                      y_vec,
                      "b-");


            std::vector<double> t(dt_vec.size());
            t[0]=0;
            std::partial_sum(dt_vec.begin(),dt_vec.end()-1,t.begin()+1);

            plt::figure(2);
            plt::named_plot("vx",t,vx_vec,"k--");
            plt::named_plot("vy",t,vy_vec,"b-");
            plt::legend();

            plt::figure(3);
            plt::named_plot("throttle",t,throttle_vec,"g-.");
            plt::named_plot("front brake",t,bf_vec,"k--");
            plt::named_plot("rear brake",t,br_vec,"b-");
            plt::named_plot("steering",t,delta_vec,"r-*");
            plt::legend();

            plt::show();
        }
        */

    private:
        std::shared_ptr<Path> path_;
        double path_width_{};
        param_t vehicle_params_,front_tire_params_,rear_tire_params_,track_params_;
        std::string database_file_ ="../output/global_planner.db";
        std::string datatable_name_ = std::string();
        std::string current_datatable_name_ = std::string();
        Dict option_;
        bool save_to_database_=false;
        DM road_constraint_upper_bound_,road_constraint_lower_bound_;
        const int N_{};

        //std::shared_ptr<valid_checker_t> valid_checker_;

        constexpr static int nx = 7;
        constexpr static int nu = 4;
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



        void save(DM& dt_array,DM& x,DM& u){
            SQLite::Database    db(database_file_, SQLite::OPEN_READWRITE|SQLite::OPEN_CREATE);
            std::cout << "SQLite database file '" << db.getFilename().c_str() << "' opened successfully\n";
            current_datatable_name_ = datatable_name_;
            if(current_datatable_name_.empty()){
                const auto now = std::chrono::system_clock::now();
                time_t rawtime;
                struct tm * timeinfo;
                char buffer[40];
                time (&rawtime);
                timeinfo = localtime(&rawtime);
                strftime(buffer,sizeof(buffer),"_%d_%H_%M_%S",timeinfo);
                current_datatable_name_ = std::string(buffer);
            }


            {
                //delete table if exist
                std::string drop_statement_string = "DROP TABLE IF EXISTS " + current_datatable_name_;
                db.exec(drop_statement_string);
                std::cout << "Save data to data table: " << current_datatable_name_ << "\n";

                //create table
                std::string table_statement_string = "create table if not exists "
                                                     +current_datatable_name_
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
                                                       "front_brake REAL,"
                                                       "rear_brake REAL"
                                                       ")";
                SQLite::Statement query(db, table_statement_string);
                try {
                    query.exec();
                } catch (SQLite::Exception &e) {
                    std::cout << "Create Table "<<current_datatable_name_<<" Error\n";
                    std::cout << e.what()<<std::endl;
                    return;
                }
                std::string query_string = "INSERT INTO " +current_datatable_name_
                                           +" (dt,tau,s,n,x,y,phi,vx,vy,omega,delta,delta_dot,throttle,front_brake,rear_brake) "
                                            "VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)";
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
                    query.bind(15,double(u(3,i)));
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
