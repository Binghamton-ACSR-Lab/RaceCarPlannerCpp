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

        UnityController(double estimated_dt = 0.1, bool save_global_data=true,bool save_local_data=true,std::string database_file = "../output/global_planner.db")
            :db_manager(database_file),estimated_dt_(estimated_dt),save_global_data_(save_global_data),save_local_data_(save_global_data),
             server_(0.1){

        }

        void run() {

            server_.run(std::bind(&UnityController::waypoints_processor, this, std::placeholders::_1),
                       std::bind(&UnityController::data_processor, this, std::placeholders::_1));
        }

        json waypoints_processor(const json &value){
            auto& wapoints_json = value.at("waypoints");
            auto& pose=value.at("pose");

            std::vector<double> waypoint_x = wapoints_json.at("ptsx");
            std::vector<double> waypoint_y = wapoints_json.at("ptsy");

            double pt_x = pose.at("x");
            double pt_y = pose.at("y");
            double phi=pose.at("psi");

            waypoints_ = DM(std::vector<std::vector<double>>{waypoint_x,waypoint_y});

            std::async(std::launch::deferred,[&waypoint_x,&waypoint_y](){
                std::ofstream outfile("../data/temp.csv");
                for(auto i=0;i<waypoint_x.size();i++)
                    outfile<<std::fixed << std::setprecision(1)<<std::setw(8)<<waypoint_x[i]<<", "<<std::setw(10)<<waypoint_y[i]<<std::endl;
                outfile.close();
            });


            path_ptr_ = std::make_shared<Path>(waypoints_);
            //std::string database_file{"../output/global_planner.db"};
            std::string vehicle_file = "../data/params/racecar.json";
            if(!std::filesystem::exists(vehicle_file)){
                std::cout<<vehicle_file<<" does not exist\n";
                return R"({"foo": "bar"})"_json;
            }
            std::ifstream ifs(vehicle_file);
            auto vehicle_params = json::parse(ifs);

            double track_width = value.at("track_width");
            int N = 200;
            GlobalPlanner optimizer(path_ptr_,track_width,vehicle_params);

            auto dm_pt = DM{pt_x,pt_y};
            auto tau0 = path_ptr_->xy2t(dm_pt);
            auto n0 =path_ptr_->f_xy_to_tn(DMVector{dm_pt,tau0})[0];

            param_t init{{"t",double(tau0)},{"n",double(n0)},{"phi",phi},{"v",0.15}};
            auto return_value = optimizer.make_plan(init,1000,N);
            if(return_value.first){
                auto& data = return_value.second;
                auto& dt = std::get<0>(data);
                auto& x = std::get<1>(data);
                auto& u = std::get<2>(data);
                global_trajectory_ptr_ = std::make_shared<GlobalTrajectory>(path_ptr_,dt,x,u);

                local_planner_ptr_ = std::make_shared<LocalPlanner>(global_trajectory_ptr_, vehicle_params,10, 0.1);

                std::future<void> save_data_thread;

                if(save_global_data_){
                    save_data_thread = std::async(std::launch::async,[&](){
                        optimizer.save(db_manager, dt,x,u);
                    });
                }


                auto ref_dm = path_ptr_->f_tn_to_xy(DMVector{x(0,Slice()),x(1,Slice())})[0];
                auto ref_x = ref_dm(0,Slice());
                auto ref_y = ref_dm(1,Slice());
                json j;
                j["ref_x"] = std::vector(ref_x->begin(),ref_x->end());
                j["ref_y"] = std::vector(ref_y->begin(),ref_y->end());
                if(save_global_data_){
                    save_data_thread.wait();
                }
                return j;
            }

            return R"({"foo": "bar"})"_json;
        }

        json data_processor(const json &value){
            //std::cout<<value<<std::endl;
            DM x0 = DM{value.at("x"),value.at("y"),value.at("psi"),value.at("vx")};
            //std::cout<<"received_x0"<<x0<<std::endl;
            auto estimate_state = x0 + estimated_dt_ * local_planner_ptr_->dynamics_model_cartesian(x0,sent_control_);
            //std::cout<<"estimated_x0"<<estimate_state<<std::endl;

            DM x,u;
            auto local_planner = local_planner_ptr_->make_plan(estimate_state);
            if(local_planner.first) {
                x = local_planner.second.first;
                u = local_planner.second.second;
                sent_control_ = u(Slice(), 0);
            }else{
                return R"({"foo": "bar"})"_json;
            }

            json j;
            auto mpc_x = x(0,Slice());
            auto mpc_y = x(1,Slice());
            j["mpc_x"] = std::vector(mpc_x->begin(),mpc_x->end());
            j["mpc_y"] = std::vector(mpc_y->begin(),mpc_y->end());
            j["throttle"] = double(u(1,0));
            j["steering_angle"] = double(u(0,0));
            return j;

        }

    private:
        WsServer server_;
        DM waypoints_;
        std::shared_ptr<Path> path_ptr_;
        std::shared_ptr<GlobalTrajectory> global_trajectory_ptr_;
        std::shared_ptr<LocalPlanner> local_planner_ptr_;
        bool save_global_data_,save_local_data_;
        //SQLite::Database    db_;
        DbManager db_manager;
        double estimated_dt_;

        DM sent_control_=DM::zeros(LocalPlanner::nu_);



    };
}

#endif //RACECARPLANNER_UNITY_CONTROLLER_H
