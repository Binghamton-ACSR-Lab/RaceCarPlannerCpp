//
// Created by acsr on 11/16/22.
//

#ifndef RACECARPLANNER_UNITY_CONTROLLER_H
#define RACECARPLANNER_UNITY_CONTROLLER_H
#include "websocket_server.hpp"
#include "local_planner/local_planner.hpp"
#include "global_planner_optimizer.hpp"
#include "db_manager.hpp"

namespace acsr {

    template<class GlobalPlanner=BicycleKinematicOptimizer,class LocalPlanner=BicycleKinematicController>
    struct UnityController {

        UnityController(bool save_global_data=true,bool save_local_data=true,std::string database_file = "../output/global_planner.db")
            :db_manager(database_file){

        }

        void run() {

            server.run(std::bind(&UnityController::waypoints_processor, this, std::placeholders::_1),
                       std::bind(&UnityController::data_processor, this, std::placeholders::_1));
        }

        void waypoints_processor(const json &value){
            std::vector<double> waypoint_x = value["waypoints"]["ptsx"];
            std::vector<double> waypoint_y = value["waypoints"]["ptsy"];

            double pt_x = value["pose"]["x"];
            double pt_y = value["pose"]["y"];
            double phi=value["pose"]["psi"];

            waypoints = DM(std::vector<std::vector<double>>{waypoint_x,waypoint_y});
            std::ofstream outfile("../data/temp.csv");


            for(auto i=0;i<waypoint_x.size();i++)
                outfile<<std::fixed << std::setprecision(1)<<std::setw(8)<<waypoint_x[i]<<", "<<std::setw(10)<<waypoint_y[i]<<std::endl;
            outfile.close();

            path_ptr_ = std::make_shared<Path>(waypoints);



            //std::string database_file{"../output/global_planner.db"};
            std::string vehicle_file = "../data/params/racecar.json";
            if(!std::filesystem::exists(vehicle_file)){
                std::cout<<vehicle_file<<" does not exist\n";
                return;
            }
            std::ifstream ifs(vehicle_file);
            auto vehicle_params = json::parse(ifs);
            /*std::string vehicle_file = "../data/params/racecar.yaml";
            if(!std::filesystem::exists(vehicle_file)){
                std::cout<<vehicle_file<<" does not exist\n";
                return;
            }

            param_t vehicle_params;
            auto vehicle_yaml = YAML::LoadFile(vehicle_file);
            std::cout<<"load vehicle config file.. Total node: "<<vehicle_yaml.size()<<std::endl;
            for(auto it = vehicle_yaml.begin();it!=vehicle_yaml.end();++it){
                vehicle_params[it->first.as<std::string>()]=it->second.as<double>();
            }*/

            double track_width = value["track_width"];
            int N = 1000;

            GlobalPlanner optimizer(path_ptr_,track_width,vehicle_params);

            auto dm_pt = DM{pt_x,pt_y};
            auto tau0 = path_ptr_->xy2t(dm_pt);
            auto n0 =path_ptr_->f_xy_to_tn(DMVector{dm_pt,tau0})[0];

            param_t init{{"t",double(tau0)},{"n",double(n0)},{"phi",phi},{"v",0.15}};
            auto return_value = optimizer.make_plan(init,-1,N);
            if(return_value.first){
                auto& data = return_value.second;
                auto& dt = std::get<0>(data);
                auto& x = std::get<1>(data);
                auto& u = std::get<2>(data);
                global_trajectory_ptr_ = std::make_shared<GlobalTrajectory>(path_ptr_,dt,x,u);

                local_planner_ptr_ = std::make_shared<LocalPlanner>(global_trajectory_ptr_, vehicle_params,10, 0.1);
                if(save_global_data_){
                    optimizer.save(db_manager, dt,x,u);//,std::vector<std::string>{"t","n","phi","v"}, std::vector<std::string>{"delta","d"});
                }
            }


        }

        json data_processor(const json &value){
            std::cout<<value<<std::endl;
            /*
            received_phi = float(telemetry['psi'])
            received_ptx = float(telemetry['x'])
            received_pty = float(telemetry['y'])
            received_vx = float(telemetry['vx'])
            received_vy = float(telemetry['vy'])*/
            return R"({"foo": "bar"})"_json;

        }





    private:
        WsServer server;
        DM waypoints;
        std::shared_ptr<Path> path_ptr_;
        std::shared_ptr<GlobalTrajectory> global_trajectory_ptr_;
        std::shared_ptr<LocalPlanner> local_planner_ptr_;
        bool save_global_data_,save_local_data_;
        //SQLite::Database    db_;
        DbManager db_manager;



    };
}

#endif //RACECARPLANNER_UNITY_CONTROLLER_H
