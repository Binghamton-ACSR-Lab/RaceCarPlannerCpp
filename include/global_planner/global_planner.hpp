//
// Created by acsr on 9/7/22.
//

#ifndef RACECARPLANNER_GLOBAL_PLANNER_HPP
#define RACECARPLANNER_GLOBAL_PLANNER_HPP

#include "track.hpp"
#include <yaml-cpp/yaml.h>
#include <filesystem>
#include "dynamics.hpp"
#include "tire_model.hpp"
#include <SQLiteCpp/SQLiteCpp.h>
#include <boost/format.hpp>
#include <filesystem>

namespace acsr{
    using param_t = std::map<std::string,double>;
    class BicycleDynamicsGlobalPlanner{
    public:
        BicycleDynamicsGlobalPlanner() = default;
        BicycleDynamicsGlobalPlanner(const std::string& track_file,const std::string& vehicle_file,const std::string& front_tire_file,const std::string& rear_tire_file){
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

            std::filesystem::path p{track_file};
            std::filesystem::path track_config_file = p.replace_extension("yaml");
            auto track_yaml = YAML::LoadFile(std::string{track_config_file});
            std::cout<<"load track config.. Total node: "<<track_yaml.size()<<std::endl;

            auto track_width = track_yaml["width"].as<double>();
            auto closed = track_yaml["closed"].as<bool>();
            track = std::make_shared<Track>(track_file,track_width,closed);
            option["max_iter"] = 600000;
            option["tol"]=1e-6;
            option["linear_solver"]="ma27";
            //option["print_time"]=true;
        }

        void plan(double start_s = 0, double end_s = -1,double n0 =0, double v0=0.15,int N = 100,bool save_to_database = true){
            auto all = Slice();
            auto _N_1 = Slice(0,N);

            auto front_tire_model = std::make_shared<PacejkaSimpleModel>(front_tire_params);
            auto rear_tire_model = std::make_shared<PacejkaSimpleModel>(rear_tire_params);

            auto cd =vehicle_params.at("Cd");
            auto cm0 = vehicle_params.at("Cm0");
            auto cm1 = vehicle_params.at("Cm1");
            auto cm2 = vehicle_params.at("Cm2");
            auto cbf = vehicle_params.at("Cbf");
            auto cbr = vehicle_params.at("Cbr");
            auto lf = vehicle_params.at("lf");
            auto lr = vehicle_params.at("lr");
            auto mass = vehicle_params.at("m");
            auto Iz = vehicle_params.at("Iz");
            auto Fz = casadi::DM::ones(2) * vehicle_params.at("m") / 2 * 9.81;

            auto s_val = casadi::DM::linspace(start_s,end_s,N+1).T();
            auto tau_array = track->s_to_t_lookup(s_val)[0];
            auto kappa_array = track->f_kappa(tau_array)[0];
            auto tangent_vec_array= track->f_tangent_vec(tau_array)[0];
            auto phi_array = DM::atan2(tangent_vec_array(1,Slice()),tangent_vec_array(0,Slice()));
            auto tangent_vec_norm = DM::sqrt((tangent_vec_array(0,Slice())*tangent_vec_array(0,Slice())+tangent_vec_array(1,Slice())*tangent_vec_array(1,Slice())));

            tau_array.to_file("/home/acsr/Documents/data/tau_array_cpp.txt");
            kappa_array.to_file("/home/acsr/Documents/data/kappa_array_cpp.txt");
            tangent_vec_array.to_file("/home/acsr/Documents/data/tangent_vec_array_cpp.txt");
            phi_array.to_file("/home/acsr/Documents/data/phi_array_cpp.txt");
            tangent_vec_norm.to_file("/home/acsr/Documents/data/tangent_vec_norm_cpp.txt");

            casadi::Opti opti;
            auto X = opti.variable(nx,N+1);
            auto U = opti.variable(nu, N);
            auto X_dot = opti.variable(nx, N);
            auto dt_sym_array = opti.variable(1,N);

            auto n_sym_array = X(IDX_X_n,_N_1);
            auto phi_sym_array = X(IDX_X_phi, _N_1);
            auto vx_sym_array = X(IDX_X_vx,_N_1);
            auto vy_sym_array = X(IDX_X_vy,_N_1);
            auto omega_sym_array = X(IDX_X_omega,_N_1);
            auto delta_sym_array = X(IDX_X_delta,_N_1);

            auto delta_dot_sym_array = U(IDX_U_Ddelta,all);
            auto throttle_sym_array = U(IDX_U_Throttle,all);
            auto front_brake_sym_array = U(IDX_U_Fb,all);
            auto rear_brake_sym_array = U(IDX_U_Rb,all);


            //auto n_sym_array = X(1,Slice());
            //auto n_obj = (casadi::MX::atan(5 * ( n_sym_array*n_sym_array - track->get_width()*track->get_width() / 4) ) + casadi::pi / 2) * 120;


            auto X0 = casadi::DM::vertcat({tau_array(0,0), n0, phi_array(0,0), v0, 0, 0, 0});
            std::cout<<"X0: "<<X0<<std::endl;

            auto dphi_c_sym_array = phi_sym_array - phi_array(0,_N_1);
            auto fx_f_sym_array = cd * throttle_sym_array - cm0 - cm1 * vx_sym_array * vx_sym_array - cm2 * vy_sym_array * vy_sym_array - cbf * front_brake_sym_array;
            auto fx_r_sym_array = cd * throttle_sym_array - cm0 - cm1 * vx_sym_array * vx_sym_array - cm2 * vy_sym_array * vy_sym_array - cbr * rear_brake_sym_array;


            auto alpha_f = -casadi::MX::atan2(omega_sym_array * lf + vy_sym_array, vx_sym_array) + delta_sym_array;
            auto alpha_r = casadi::MX::atan2(omega_sym_array * lr - vy_sym_array, vx_sym_array);

            auto fy_f_sym_array = front_tire_model->get_lateral_force(alpha_f, double(Fz(0)));
            auto fy_r_sym_array = rear_tire_model->get_lateral_force(alpha_r, double(Fz(1)));

            X_dot(0, all) = dt_sym_array*(vx_sym_array * MX::cos(dphi_c_sym_array) - vy_sym_array * MX::sin(dphi_c_sym_array))/(tangent_vec_norm(0,_N_1)*(1-n_sym_array*kappa_array(0,_N_1)));
            X_dot(1, all) = dt_sym_array*(vx_sym_array * MX::sin(dphi_c_sym_array) + vy_sym_array* MX::cos(dphi_c_sym_array)) ;
            X_dot(2, all) = dt_sym_array*omega_sym_array;
            X_dot(3, all) = dt_sym_array * (fx_r_sym_array + fx_f_sym_array * MX::cos(delta_sym_array) - fy_f_sym_array * MX::sin(delta_sym_array) + mass * vy_sym_array* omega_sym_array)/ mass;
            X_dot(4, all) = dt_sym_array * (fy_r_sym_array + fx_f_sym_array * MX::sin(delta_sym_array) + fy_f_sym_array * MX::cos(delta_sym_array) - mass * vx_sym_array * omega_sym_array)/ mass;
            X_dot(5, all) = dt_sym_array * (fy_f_sym_array * lf * MX::cos(delta_sym_array) + fx_f_sym_array * lf * MX::sin(delta_sym_array) - fy_r_sym_array * lr)/ Iz ;
            X_dot(6, all) = dt_sym_array * delta_dot_sym_array;

            auto track_width = track->get_width();
            auto n_obj = (MX::atan(5 * (n_sym_array * n_sym_array - (track_width / 2) *(track_width / 2))) + casadi::pi / 2) * 12;
            opti.minimize(MX::sum2(dt_sym_array) + MX::dot(n_obj, n_obj));
            //opti.minimize(MX::sum2(dt_sym_array));
            //inital conditions
            opti.subject_to(X(all, 0) == X0);
            opti.subject_to(X(0, all) == tau_array);
            opti.subject_to(dt_sym_array >0);
            //dynamics
            opti.subject_to(X(all, Slice(1,N+1)) == X(all, _N_1) + X_dot);

            /*
            opti.subject_to(X(0,Slice(1,N+1))==X(0,Slice(0,N))+dt_sym_array*X_dot(0,Slice()));
            opti.subject_to(X(1,Slice(1,N+1))==X(1,Slice(0,N))+dt_sym_array*X_dot(1,Slice()));
            opti.subject_to(X(2,Slice(1,N+1))==X(2,Slice(0,N))+dt_sym_array*X_dot(2,Slice()));
            opti.subject_to(X(3,Slice(1,N+1))==X(3,Slice(0,N))+dt_sym_array*X_dot(3,Slice()));
            opti.subject_to(X(4,Slice(1,N+1))==X(4,Slice(0,N))+dt_sym_array*X_dot(4,Slice()));
            opti.subject_to(X(5,Slice(1,N+1))==X(5,Slice(0,N))+dt_sym_array*X_dot(5,Slice()));
            opti.subject_to(X(6,Slice(1,N+1))==X(6,Slice(0,N))+dt_sym_array*X_dot(6,Slice()));
            */


            //state boundary
            opti.subject_to(opti.bounded(vehicle_params.at("delta_min"), X(IDX_X_delta,all), vehicle_params.at("delta_max")));

            //control boundary
            opti.subject_to(opti.bounded(vehicle_params.at("delta_dot_min"), delta_dot_sym_array, vehicle_params.at("delta_dot_max")));
            opti.subject_to(opti.bounded(0, throttle_sym_array, 1));
            opti.subject_to(opti.bounded(0, front_brake_sym_array, 1));
            opti.subject_to(opti.bounded(0, rear_brake_sym_array, 1));

            //opti.subject_to(opti.bounded(-track->get_width()/2,X(IDX_X_n,all),track->get_width()/2));

            auto X_guess = casadi::DM::zeros(nx,N+1);
            X_guess(0,Slice()) = tau_array;
            X_guess(2,Slice()) = phi_array;
            X_guess(3,Slice()) = v0;
            //for(auto i=0;i<N+1;++i){
            //    X_guess(Slice(), i) = X0;
            //}
            opti.set_initial(X, X_guess);

            opti.solver("ipopt", Dict(), option);
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


        void make_plan(double start_s = 0, double end_s = -1,double n0 =0, double v0=0.15,int N = 100,bool save_to_database = true){
            auto front_tire_model = std::make_shared<PacejkaSimpleModel>(front_tire_params);
            auto rear_tire_mode = std::make_shared<PacejkaSimpleModel>(rear_tire_params);
            BicycleDynamicsByParametricArc<PacejkaSimpleModel,PacejkaSimpleModel> dynamics(vehicle_params,track,front_tire_model,rear_tire_mode);

            auto s_val = casadi::DM::linspace(start_s,end_s,N+1).T();
            auto tau_array = track->s_to_t_lookup(s_val)[0];

            casadi::Opti opti;
            auto X = opti.variable(dynamics.nx(),N+1);
            //auto X_dot = opti.variable(dynamics.nx(), N);
            auto U = opti.variable(dynamics.nu(), N);
            auto dt_sym_array = opti.variable(1,N);

            auto X_dot = dynamics.updata(X,U);
            std::cout<<"x_dot size: "<<X_dot.size()<<std::endl;

            auto n_sym_array = X(1,Slice());
            //auto n_obj = (casadi::MX::atan(5 * ( n_sym_array*n_sym_array - track->get_width()*track->get_width() / 4) ) + casadi::pi / 2) * 120;
            opti.minimize(casadi::MX::sum2(dt_sym_array));

            //auto tau0 = track->s_to_t_lookup(DM{start_s});
            auto phi_array = track->f_phi(tau_array)[0];
            auto X0 = casadi::DM::vertcat({tau_array(0,0), n0, phi_array(0,0), v0, 0, 0, 0});
            std::cout<<"X0: "<<X0<<std::endl;
            //inital conditions
            opti.subject_to(X(Slice(), 0) == X0);
            opti.subject_to(X(0, Slice()) == tau_array);
            //dynamics
            opti.subject_to(X(0,Slice(1,N+1))==X(0,Slice(0,N))+dt_sym_array*X_dot(0,Slice()));
            opti.subject_to(X(1,Slice(1,N+1))==X(1,Slice(0,N))+dt_sym_array*X_dot(1,Slice()));
            opti.subject_to(X(2,Slice(1,N+1))==X(2,Slice(0,N))+dt_sym_array*X_dot(2,Slice()));
            opti.subject_to(X(3,Slice(1,N+1))==X(3,Slice(0,N))+dt_sym_array*X_dot(3,Slice()));
            opti.subject_to(X(4,Slice(1,N+1))==X(4,Slice(0,N))+dt_sym_array*X_dot(4,Slice()));
            opti.subject_to(X(5,Slice(1,N+1))==X(5,Slice(0,N))+dt_sym_array*X_dot(5,Slice()));
            opti.subject_to(X(6,Slice(1,N+1))==X(6,Slice(0,N))+dt_sym_array*X_dot(6,Slice()));

            opti.subject_to(dt_sym_array >0);

            //state boundary
            opti.subject_to(opti.bounded(vehicle_params.at("delta_min"), X(6,Slice()), vehicle_params.at("delta_max")));
            opti.subject_to(opti.bounded(-track->get_width()/2,n_sym_array,track->get_width()/2));
            //control boundary
            opti.subject_to(opti.bounded(vehicle_params.at("delta_dot_min"), U(0, Slice()), vehicle_params.at("delta_dot_max")));
            opti.subject_to(opti.bounded(0, U(1, Slice()), 1));
            opti.subject_to(opti.bounded(0, U(2, Slice()), 1));
            opti.subject_to(opti.bounded(0, U(3, Slice()), 1));

            auto X_guess = casadi::DM::zeros(dynamics.nx(),N+1);
            //X_guess(0,Slice()) = tau_array;
            //X_guess(2,Slice()) = phi_array;
            //X_guess(3,Slice()) = v0;
            for(auto i=0;i<N+1;++i){
                X_guess(Slice(), i) = X0;
            }
            opti.set_initial(X, X_guess);

            opti.solver("ipopt", Dict(), option);
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
        std::shared_ptr<Track> track;
        param_t vehicle_params,front_tire_params,rear_tire_params,track_params;
        std::string database_file ="../output/global_planner.db";
        std::string datatable_name = std::string();
        Dict option;

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


            SQLite::Database    db(database_file);
            std::string table = datatable_name;

            if(table.empty()){
                const auto now = std::chrono::system_clock::now();
                time_t rawtime;
                struct tm * timeinfo;
                char buffer[40];
                time (&rawtime);
                timeinfo = localtime(&rawtime);
                strftime(buffer,sizeof(buffer),"_%d_%H_%M_%S",timeinfo);
                table = std::string(buffer);
            }

            std::cout << "SQLite database file '" << db.getFilename().c_str() << "' opened successfully\n";
            {
                //delete table if exist
                std::string drop_statement_string = "DROP TABLE IF EXISTS " + table;
                db.exec(drop_statement_string);

                //create table
                std::string table_statement_string = "create table if not exists "
                                                     +table
                                                     + "("
                                                       "id INTEGER PRIMARY KEY AUTOINCREMENT, "
                                                       "dt REAL NOT NULL,"
                                                       "tau REAL NOT NULL,"
                                                       "s REAL NOT NULL,"
                                                       "n REAL NOT NULL,"
                                                       "x REAL NOT NULL,"
                                                       "y REAL NOT NULL,"
                                                       "phi REAL NOT NULL,"
                                                       "vx REAL NOT NULL,"
                                                       "vy REAL NOT NULL,"
                                                       "omega REAL NOT NULL,"
                                                       "delta REAL NOT NULL,"
                                                       "delta_dot REAL,"
                                                       "throttle REAL,"
                                                       "front_brake REAL,"
                                                       "rear_brake REAL"
                                                       ")";
                SQLite::Statement query(db, table_statement_string);
                try {
                    query.exec();
                } catch (SQLite::Exception &e) {
                    std::cout << "Create Table "<<table<<" Error\n";
                    std::cout << e.what()<<std::endl;
                    return;
                }

                std::string query_string = "INSERT INTO " +table
                        +" (dt,tau,s,n,x,y,phi,vx,vy,omega,delta,delta_dot,throttle,front_brake,rear_brake) "
                         "VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)";
                query = SQLite::Statement(db, query_string);

                auto N = x.columns();
                auto dm_xy = track->f_sn_to_xy(std::vector<DM>{x(0),x(1)})[0];
                auto s_array = track->t_to_s_lookup(x(0))[0];
                for(auto i=0;i<N;++i){
                    query.bind(1,double(dt_array(0,i)));
                    query.bind(2,double(s_array(0,i)));
                    query.bind(3,double(x(1,i)));
                    query.bind(4,double(x(1,i)));
                    query.bind(5,double(dm_xy(0,i)));
                    query.bind(6,double(dm_xy(1,i)));
                    query.bind(7,double(x(2,i)));
                    query.bind(8,double(x(3,i)));
                    query.bind(9,double(x(4,i)));
                    query.bind(10,double(x(5,i)));
                    query.bind(11,double(x(6,i)));
                    query.bind(12,double(u(1,i)));
                    query.bind(13,double(u(2,i)));
                    query.bind(14,double(u(3,i)));
                    query.bind(15,double(u(4,i)));
                    try {
                        query.exec();
                    } catch (SQLite::Exception &e) {
                        std::cout << "Insert Solution Error\n";
                        std::cout << e.what();
                        return;
                    }
                }
            }
        }

    };
}

#endif //RACECARPLANNER_GLOBAL_PLANNER_HPP
