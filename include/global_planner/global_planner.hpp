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

namespace acsr{
    using param_t = std::map<std::string,double>;
    class BicycleDynamicsGlobalPlanner{
    public:
        BicycleDynamicsGlobalPlanner() = default;
        BicycleDynamicsGlobalPlanner(const std::string& track_file,const std::string& vehicle_file,const std::string& front_tire_file,const std::string& rear_tire_file){

            //load params
            auto vehicle_yaml = YAML::Load(vehicle_file);
            for(auto it = vehicle_yaml.begin();it!=vehicle_yaml.end();++it){
                vehicle_params[it->first.as<std::string>()]=it->second.as<double>();
            }

            auto front_tire_yaml = YAML::Load(front_tire_file);
            for(auto it = front_tire_yaml.begin();it!=front_tire_yaml.end();++it){
                front_tire_params[it->first.as<std::string>()]=it->second.as<double>();
            }

            auto rear_tire_yaml = YAML::Load(rear_tire_file);
            for(auto it = rear_tire_yaml.begin();it!=rear_tire_yaml.end();++it){
                rear_tire_params[it->first.as<std::string>()]=it->second.as<double>();
            }

            std::filesystem::path p{track_file};
            std::filesystem::path track_config_file = p.replace_extension("yaml");
            auto track_yaml = YAML::Load(std::string{track_config_file});


            track = std::make_shared<Track>(track_file,track_yaml["width"].as<double>(),track_yaml["closed"].as<bool>());
        }

        void make_plan(double start_s = 0,double end_s = -1,int resolution = 100){
            auto front_tire_model = std::make_shared<PacejkaSimpleModel>(front_tire_params);
            auto rear_tire_mode = std::make_shared<PacejkaSimpleModel>(rear_tire_params);
            BicycleDynamicsByParametricArc<PacejkaSimpleModel,PacejkaSimpleModel> dynamics(vehicle_params,track,front_tire_model,rear_tire_mode);

            casadi::Opti opti;
            auto X = opti.variable();

        }


    private:
        std::shared_ptr<Track> track;
        param_t vehicle_params,front_tire_params,rear_tire_params,track_params;
    };
}

#endif //RACECARPLANNER_GLOBAL_PLANNER_HPP
