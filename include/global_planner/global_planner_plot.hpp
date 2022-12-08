//
// Created by acsr on 12/1/22.
//

#ifndef RACECARPLANNER_GLOBAL_PLANNER_PLOT_HPP
#define RACECARPLANNER_GLOBAL_PLANNER_PLOT_HPP
#include "global_trajectory.hpp"
#include "csv_reader.hpp"
#include "SQLiteCpp/SQLiteCpp.h"
#include "matplotlibcpp.h"

namespace acsr {

    template<int nx,int nu,int vx_index,int vy_index=-1>
    struct GlobalPlannerPlotter {
        GlobalPlannerPlotter(std::shared_ptr<GlobalTrajectory> trajectory_ptr):trajectory_ptr_(trajectory_ptr) {
        }

        GlobalPlannerPlotter(std::shared_ptr<Path> path_ptr, const std::string& trajectory_table=std::string{},const std::string& database_file="../output/global_planner.db"){

            SQLite::Database db(database_file);
            SQLite::Statement query1(db,"select name from sqlite_sequence");

            std::vector<double> dt_vec;
            std::vector<std::vector<double>> x_vec,u_vec;


            std::string table_name;
            while (query1.executeStep()){
                table_name = query1.getColumn(0).getString();
            }

            SQLite::Statement query2(db,"select * from "+table_name);
            while (query2.executeStep()){
                dt_vec.push_back(query2.getColumn(0));
                std::vector<double> x,u;
                for(auto i=0;i<nx;++i){
                    x.push_back(query2.getColumn(1+i));
                }
                for(auto i=0;i<nu;++i){
                    u.push_back(query2.getColumn(1+nx+i));
                }
                x_vec.push_back(x);
            }
            trajectory_ptr_ = std::make_shared<GlobalTrajectory>(path_ptr,DM(dt_vec),DM(x_vec).T(),DM(u_vec).T());
        }

        GlobalPlannerPlotter(const std::string& track_file="../data/temp.csv",const std::string& trajectory_table=std::string{},const std::string& database_file="../output/global_planner.db"){
            DM waypoints;
            CSVReader reader(track_file);
            auto csv_data = reader.read();
            if(!csv_data.toDM(waypoints)){
                std::cout<<"Convert csv data to dm fail\n";
                return;
            }
            auto path_ptr = std::make_shared<Path>(waypoints);
            SQLite::Database db(database_file);

            std::vector<double> dt_vec;
            std::vector<std::vector<double>> x_vec,u_vec;


            std::string table_name = trajectory_table;
            if(table_name.empty()) {
                SQLite::Statement query1(db, "select name from sqlite_sequence");
                while (query1.executeStep()) {
                    table_name = query1.getColumn(0).getString();
                }
            }
            dt_vec.push_back(0);
            SQLite::Statement query2(db,"select * from "+table_name);
            while (query2.executeStep()){
                dt_vec.push_back(query2.getColumn(1));
                std::vector<double> x,u;
                for(auto i=0;i<nx;++i){
                    x.push_back(query2.getColumn(2+i));
                }
                for(auto i=0;i<nu;++i){
                    u.push_back(query2.getColumn(2+nx+i));
                }
                x_vec.push_back(x);
                u_vec.push_back(u);
            }
            dt_vec.pop_back();
            u_vec.pop_back();
            trajectory_ptr_ = std::make_shared<GlobalTrajectory>(path_ptr,DM(dt_vec),DM(x_vec).T(),DM(u_vec).T());
        }

        void plot(double width, int N = 500){
            namespace plt = matplotlibcpp;

            auto path_ptr = trajectory_ptr_->get_path();
            auto t_max = path_ptr->get_max_tau();
            auto ts = casadi::DM::linspace(0, t_max, N * t_max);
            DM ns_zeros = DM::zeros(ts.rows(),ts.columns());
            DM ns_ones = DM::ones(ts.rows(),ts.columns());
            auto center_line = path_ptr->f_tn_to_xy(DMVector {ts.T(),ns_zeros.T()})[0];
            auto inner_line = path_ptr->f_tn_to_xy(std::vector<DM>{ts.T(),+width/2*ns_ones.T()})[0];
            auto outer_line = path_ptr->f_tn_to_xy(std::vector<DM>{ts.T(),-width/2*ns_ones.T()})[0];
            auto center_x_dm = center_line(0,Slice());
            auto center_y_dm = center_line(1,Slice());
            auto inner_x_dm = inner_line(0,Slice());
            auto inner_y_dm = inner_line(1,Slice());
            auto outer_x_dm = outer_line(0,Slice());
            auto outer_y_dm = outer_line(1,Slice());

            auto dm_t = trajectory_ptr_->get_t_vec().T();
            auto state = trajectory_ptr_->f_t_to_state(DMVector{dm_t})[0];
            auto trajectory_line = path_ptr->f_tn_to_xy(DMVector {state(0,Slice()),state(1,Slice())})[0];
            auto trajectory_x_dm = trajectory_line(0,Slice());
            auto trajectory_y_dm = trajectory_line(1,Slice());

            plt::plot(std::vector<double>{center_x_dm->begin(),center_x_dm->end()},
                      std::vector<double>{center_y_dm->begin(),center_y_dm->end()},
                      "r-",
                      std::vector<double>{inner_x_dm->begin(),inner_x_dm->end()},
                      std::vector<double>{inner_y_dm->begin(),inner_y_dm->end()},
                      "g-",
                      std::vector<double>{outer_x_dm->begin(),outer_x_dm->end()},
                      std::vector<double>{outer_y_dm->begin(),outer_y_dm->end()},
                      "b-",
                      std::vector<double>{trajectory_x_dm->begin(),trajectory_x_dm->end()},
                      std::vector<double>{trajectory_y_dm->begin(),trajectory_y_dm->end()},
                      "k--");

            plt::show();
        }

    private:
        std::shared_ptr<GlobalTrajectory> trajectory_ptr_;

    };


}


#endif //RACECARPLANNER_GLOBAL_PLANNER_PLOT_H
