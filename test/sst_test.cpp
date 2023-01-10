//
// Created by xli185 on 12/28/22.
//
#include "test_map.hpp"
#include "sst.hpp"
#include "dynamics.hpp"


using namespace acsr;
using nlohmann::json;

int main(){
    //read planner params files
    std::string planner_file = "../data/params/planner.json";
    if(!std::filesystem::exists(planner_file)){
        std::cout<<planner_file<<" does not exist\n";
        return 1;
    }
    std::ifstream ifs_planner(planner_file);
    auto planner_params = json::parse(ifs_planner);


    //read racecar params files
    std::string racecar_file = "../data/params/sst_test_racecar.json";
    if(!std::filesystem::exists(racecar_file)){
        std::cout<<racecar_file<<" does not exist\n";
        return 1;
    }
    std::ifstream ifs_racecar(planner_file);
    auto racecar_params = json::parse(ifs_racecar);
    //initialize dynamics model
    auto dynamics = std::make_shared<BicycleKinematicNoAcc>(racecar_params);

    auto test_map = std::make_shared<TestMap>("../data/map/test_map.xml",0.02);

    SST<BicycleKinematicNoAcc,TestMap> sst(dynamics,test_map,0.5,planner_params);

    return 0;
}