#include <iostream>

#include <casadi/casadi.hpp>
#include "include/utility/utility.hpp"
#include <gtest/gtest.h>
#include "track.hpp"
#include "global_planner.hpp"
#include <yaml-cpp/yaml.h>
#include "matplotlibcpp.h"

using namespace acsr;
template<class T>
T update(T& x,T& u){
    return T::vertcat({x,u});
}

void planner_test(){
    std::string track_data_file = "../data/tracks/temp_nwh.csv";
    std::string track_config_file = "../data/params/racecar_nwh.yaml";
    std::string front_tire_file = "../data/params/nwh_tire.yaml";
    std::string rear_tire_file = "../data/params/nwh_tire.yaml";
    BicycleDynamicsGlobalPlanner planner(track_data_file,track_config_file,front_tire_file,rear_tire_file);
    planner.plan_refine(110.01,610.01,0,0.15,200);
    planner.plot_trajectory();
}

void track_test(){
    namespace plt = matplotlibcpp;
    std::string full_track_data_file = "../data/tracks/temp_nwh.csv";
    std::string half_track_data_file = "../data/tracks/temp_nwh_half.csv";

    Track full_track(full_track_data_file,7.0,true);
    Track half_track(half_track_data_file,7.0, false);

    auto full_tau = (DM::linspace(0,full_track.get_max_tau(),100*full_track.get_max_tau()+1)).T();
    auto half_tau = (DM::linspace(0,half_track.get_max_tau(),100*half_track.get_max_tau()+1)).T();
    auto full_zeros = DM::zeros(1,full_tau.columns());
    auto half_zeros = DM::zeros(1,half_tau.columns());

    auto full_centerline = full_track.f_tn_to_xy(DMVector{full_tau,full_zeros})[0];
    auto half_centerline = half_track.f_tn_to_xy(DMVector{half_tau,half_zeros})[0];

    auto full_x = full_centerline(0,Slice());
    auto full_y = full_centerline(1,Slice());

    auto half_x = half_centerline(0,Slice());
    auto half_y = half_centerline(1,Slice());

    plt::plot(std::vector<double>(full_x->begin(),full_x->end()),
              std::vector<double>(full_y->begin(),full_y->end()),
              "k-",
              std::vector<double>(half_x->begin(),half_x->end()),
              std::vector<double>(half_y->begin(),half_y->end()),
              "r--");
    plt::show();

}

int main(int argc, char **argv) {

    track_test();




    //::testing::InitGoogleTest(&argc, argv);
    //return RUN_ALL_TESTS();
    return 0;

}
