#include <iostream>

#include <casadi/casadi.hpp>
#include "include/utility/utility.hpp"
#include <gtest/gtest.h>
#include "track.hpp"
#include "global_planner.hpp"
using namespace acsr;


template<class T>
T update(T& x,T& u){
    return T::vertcat({x,u});
}

int main(int argc, char **argv) {

    //Track track("../data/tracks/temp.csv",5,true);
    //track.plot();
    std::string track_file = "../data/tracks/temp_nwh.csv";
    std::string front_tire_file = "../data/params/nwh_tire.yaml";
    std::string& rear_tire_file = "../data/params/nwh_tire.yaml";
    BicycleDynamicsGlobalPlanner planner("")



    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();


}
