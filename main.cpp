#include <iostream>

#include <casadi/casadi.hpp>
#include "include/utility/utility.hpp"
#include <gtest/gtest.h>
#include "track.hpp"
#include "global_planner.hpp"
using namespace acsr;
#include <yaml-cpp/yaml.h>

template<class T>
T update(T& x,T& u){
    return T::vertcat({x,u});
}

int main(int argc, char **argv) {

    DM a{1,2,3,4,5};
    a = a.T();
    auto b = a(0,Slice(1,5));
    std::cout<<"b: "<<b<<std::endl;

    std::string track_data_file = "../data/tracks/temp_nwh.csv";
    std::string track_config_file = "../data/params/racecar_nwh.yaml";
    std::string front_tire_file = "../data/params/nwh_tire.yaml";
    std::string rear_tire_file = "../data/params/nwh_tire.yaml";
    BicycleDynamicsGlobalPlanner planner(track_data_file,track_config_file,front_tire_file,rear_tire_file);
    planner.plan(0,200);


    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();


}
