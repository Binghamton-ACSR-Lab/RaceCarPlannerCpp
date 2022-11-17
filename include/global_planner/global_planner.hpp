//
// Created by acsr on 11/16/22.
//

#ifndef RACECARPLANNER_GLOBAL_PLANNER_HPP
#define RACECARPLANNER_GLOBAL_PLANNER_HPP

#include "path.hpp"
#include <yaml-cpp/yaml.h>
#include <filesystem>
#include "dynamics.hpp"
#include "tire_model.hpp"
#include <SQLiteCpp/SQLiteCpp.h>
#include <boost/format.hpp>
#include <filesystem>
#include <execution>

using casadi::DM;
using casadi::MX;
namespace acsr{


    template <int nx_,int nu_,int N_>
    class AcsrGlobalPlanner{
    public:
        AcsrGlobalPlanner()=default;

        template<class DYNAMICS, class GOAL>
        std::pair<bool,std::tuple<DM,DM,DM>> make_plan(DYNAMICS update,GOAL goal,DM& x0, DM& tau_array, DM& state_boundary, DM& control_boundary,casadi::Dict solver_option={},casadi::Dict casadi_option={}){
            auto all = casadi::Slice();
            auto _0_N = casadi::Slice(0,N_);
            auto _N_1 = casadi::Slice(1,N_+1);

#ifdef DEBUG
            assert(tau_array.rows()==1);
            assert(state_boundary.rows()==2 && state_boundary.columns()==nx_);
            assert(control_boundary.rows()==2 && control_boundary.columns()==nu_);
            assert(x0.is_vector() && x0->size()==nx_);
#endif

            casadi::Opti opti;
            auto X = opti.variable(nx_,N_+1);
            auto U = opti.variable(nu_, N_);
            auto dt_sym_array = opti.variable(1,N_);

            auto k1 = update(X(all,_0_N),U(all,all));
            auto k2 = update(X(all,_0_N)+dt_sym_array/2*k1, U(all,k));
            auto k3 = update(X(all,_0_N)+dt_sym_array/2*k2, U(all,k));
            auto k4 = update(X(all,_0_N)+dt_sym_array*k3,   U(all,k));
            auto x_next = X(all,_0_N) + dt_sym_array/6*(k1+2*k2+2*k3+k4);


            opti.minimize(goal(dt_sym_array,X,U));

            //dynamics
            opti.subject_to(X(all,_0_N) + dt_sym_array/6*(k1+2*k2+2*k3+k4) == X(all,_N_1))

            //reference
            opti.subject_to(X(0, all) == tau_array);

            //initial
            opti.subject_to(X(all, 0) == X0);
            opti.subject_to(dt_sym_array >0);

            //state boundary
            for(auto i=0;i<nx_;++i){
                if(state_boundary(1,i)>state_boundary(0,i))
                    opti.subject_to(opti.bounded(state_boundary(0,i), X(i,all), state_boundary(1,i)));
            }

            //control boundary
            for(auto i=0;i<nu_;++i){
                if(control_boundary(1,i)>control_boundary(0,i))
                    opti.subject_to(opti.bounded(control_boundary(0,i), X(i,all), control_boundary(1,i)));
            }

            auto X_guess = casadi::DM::zeros(nx_,N_+1);
            for(auto i=0;i<N_+1;++i){
                X_guess(Slice(), i) = X0;
            }
            opti.set_initial(X, X_guess);

            opti.solver("ipopt", casadi_option, solver_option);
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


        std::string save(SQLite::Database &db, DM& dt_array,DM& x,DM& u,const std::vector<std::string>& x_headers={},const std::vector<std::string>& u_headers={}){
            std::cout << "SQLite database file '" << db.getFilename().c_str() << "' opened successfully\n";
            const auto now = std::chrono::system_clock::now();
            time_t rawtime;
            struct tm * timeinfo;
            char buffer[40];
            time (&rawtime);
            timeinfo = localtime(&rawtime);
            strftime(buffer,sizeof(buffer),"_%d_%H_%M_%S",timeinfo);
            auto datatable_name = std::string(buffer);

            std::ostringstream os;
            os<<"create table if not exists "<<datatable_name<<"(id INTEGER PRIMARY KEY AUTOINCREMENT, dt REAL ";
            if(x_headers.size()<nx_){
                for (int i = 0; i < nx_; ++i) {
                    os<<",x_"<<i<<" REAL";
                }
            }else{
                for (int i = 0; i < nx_; ++i) {
                    os<<","<<x_headers[i]<<" REAL";
                }
            }
            if(u_headers.size()<nu_){
                for (int i = 0; i < nu_; ++i) {
                    os<<",u_"<<i<<" REAL";
                }
            }else{
                for (int i = 0; i < nx_; ++i) {
                    os<<","<<u_headers[i]<<" REAL";
                }
            }
            os<<")";
            auto table_statement_string = os.str();

            {
                //delete table if exist
                std::string drop_statement_string = "DROP TABLE IF EXISTS " + datatable_name;
                db.exec(drop_statement_string);
                std::cout << "Save data to data table: " << datatable_name << "\n";

                SQLite::Statement query(db, table_statement_string);
                try {
                    query.exec();
                } catch (SQLite::Exception &e) {
                    std::cout << "Create Table "<<datatable_name<<" Error\n";
                    std::cout << e.what()<<std::endl;
                    return std::string{};
                }

                //create query string "INSERT INTO table VALUES(?,?,?,?)"
                os.clear();
                os << "INSERT INTO "<<datatable_name<<" VALUES(?";
                for(int i = 0; i < nx_+nu_; i++)
                    os << ",?";
                os<<")";
                auto query_string = os.str();

                //insert data rows
                for(auto i=0;i<N_;++i){
                    query = SQLite::Statement(db, query_string);
                    query.bind(1,double(dt_array(0,i)));
                    for(auto j=0;j<nx_;++j){
                        query.bind(j+2,double(x(j,i)));
                    }
                    for(auto j=0;j<nu_;++j){
                        query.bind(j+2+nx_,double(u(0,i)));
                    }

                    try {
                        query.exec();
                    } catch (SQLite::Exception &e) {
                        std::cout << "Insert Solution Error\n";
                        std::cout << e.what();
                        return std::string{};
                    }
                }

                //insert last data row
                query = SQLite::Statement(db, query_string);
                for(auto j=0;j<nx_;++j){
                    query.bind(j+2,double(x(j,N_)));
                }
                try {
                    query.exec();
                } catch (SQLite::Exception &e) {
                    std::cout << "Insert last dataset Error\n";
                    std::cout << e.what();
                    return std::string{};
                }
            }
            std::cout<<"save to database finished\n";
            return datatable_name;
        }

    private:
        //std::string datatable_name_ = std::string();
        //std::string current_datatable_name_ = std::string();



    };


}

#endif //RACECARPLANNER_GLOBAL_PLANNER_HPP
