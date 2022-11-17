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

namespace acsr{

    template <int nx_,int nu_,int N_>
    class AcsrGlobalPlanner{
    public:
        AcsrGlobalPlanner(){}

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
                for(auto i=0;i<N_-1;++i){
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
