#include <iostream>

#include <casadi/casadi.hpp>
#include "include/utility/utility.hpp"
#include <gtest/gtest.h>
#include "track.hpp"
#include "global_planner_track.hpp"

#include <yaml-cpp/yaml.h>
#include <matplotlibcpp.h>

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
    BicycleDynamicsTwoBrakeGlobalPlanner planner(track_data_file,track_config_file,front_tire_file,rear_tire_file);
    planner.plan_refine(210.01,570.01,0,0.5,120);
    planner.plot_trajectory();
}

void track_test(){
    namespace plt = matplotlibcpp;
    auto full_track_file = "../data/tracks/temp_nwh.csv";
    auto half_track_file = "../data/tracks/temp_nwh_half.csv";

    Track full_track(full_track_file,7,true);
    Track half_track(half_track_file,7,false);

    const int resolution = 100;
    auto full_tau_max = full_track.get_max_tau();
    auto full_tau = DM::linspace(0,full_tau_max,full_tau_max*resolution+1).T();
    auto full_zeros = DM::zeros(1,full_tau.columns());
    auto full_centerline = full_track.f_tn_to_xy(DMVector{full_tau,full_zeros})[0];

    auto half_tau_max = half_track.get_max_tau();
    auto half_tau = DM::linspace(0,half_tau_max,half_tau_max*resolution+1).T();
    auto half_zeros = DM::zeros(1,half_tau.columns());
    auto half_centerline = half_track.f_tn_to_xy(DMVector{half_tau,half_zeros})[0];

    auto full_center_x = full_centerline(0,Slice());
    auto full_center_y = full_centerline(1,Slice());
    auto half_center_x = half_centerline(0,Slice());
    auto half_center_y = half_centerline(1,Slice());


    plt::plot(std::vector<double>{full_center_x->begin(),full_center_x->end()},
              std::vector<double>{full_center_y->begin(),full_center_y->end()},
              "y--",
              std::vector<double>{half_center_x->begin(),half_center_x->end()},
              std::vector<double>{half_center_y->begin(),half_center_y->end()},
              "k-");
    plt::show();


}

int main(int argc, char **argv) {

    //planner_test();
    track_test();


    //::testing::InitGoogleTest(&argc, argv);
    //return RUN_ALL_TESTS();
    return 0;

}
