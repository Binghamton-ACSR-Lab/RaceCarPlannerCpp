#include <iostream>

#include <casadi/casadi.hpp>
#include "include/utility/utility.hpp"
#include <gtest/gtest.h>
#include "track.hpp"
#include "global_planner_track.hpp"
#include <yaml-cpp/yaml.h>
#include "matplotlibcpp.h"
#include <tbb/tbb.h>
#include <execution>
#include "path.hpp"
#include "test_map.hpp"
#include "path_preprocessor.hpp"

#include "test.h"

using namespace acsr;
namespace plt = matplotlibcpp;

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
    planner.plan_refine(110.01,410.01,0,0.15,100);
    planner.plot_trajectory();
}

void track_test(){
    namespace plt = matplotlibcpp;
    std::string full_track_data_file = "../data/tracks/temp_nwh.csv";
    std::string half_track_data_file = "../data/tracks/temp_nwh_half.csv";

    Track full_track(full_track_data_file,7.0,true);
    Track half_track(half_track_data_file,7.0, false);

    std::cout<<half_track.get_max_tau()<<std::endl;
    auto full_tau = (DM::linspace(0,full_track.get_max_tau()-1e-6,100*full_track.get_max_tau()+1)).T();
    auto half_tau = (DM::linspace(0,half_track.get_max_tau()-1e-6,100*half_track.get_max_tau()+1)).T();
    auto full_zeros = DM::zeros(1,full_tau.columns());
    auto half_zeros = DM::zeros(1,half_tau.columns());

    auto full_centerline = full_track.f_tn_to_xy(DMVector{full_tau,full_zeros})[0];
    auto half_centerline = half_track.f_tn_to_xy(DMVector{half_tau,half_zeros})[0];

    auto full_x = full_centerline(0,Slice());
    auto full_y = full_centerline(1,Slice());

    auto half_x = half_centerline(0,Slice());
    auto half_y = half_centerline(1,Slice());

    std::cout<<half_centerline(Slice(),Slice(0,100))<<std::endl;
    std::cout<<"------------\n";
    std::cout<<half_centerline(Slice(),Slice(half_centerline.columns()-110,half_centerline.columns()))<<std::endl;



    plt::plot(std::vector<double>(full_x->begin(),full_x->end()),
              std::vector<double>(full_y->begin(),full_y->end()),
              "k-",
              std::vector<double>(half_x->begin(),half_x->end()),
              std::vector<double>(half_y->begin(),half_y->end()),
              "r--");
    plt::show();

}

void path_preprocessor_test(){
    std::vector<std::vector<double>> pts;
    auto filename = std::string("../data/map/test_map.xml");
    tinyxml2::XMLDocument doc;
    doc.LoadFile( filename.c_str());
    auto root = doc.RootElement();
    for (tinyxml2::XMLElement* element = root->FirstChildElement(); element != nullptr; element = element->NextSiblingElement())
    {
        if(strcmp(element->Name(),"waypoints")==0){
            for (tinyxml2::XMLElement* child = element->FirstChildElement(); child != nullptr; child = child->NextSiblingElement()){
                if(strcmp(child->Name(),"point")==0) {
                    auto x = child->DoubleAttribute("x");
                    auto y = child->DoubleAttribute("y");
                    pts.push_back({x, y});
                }
            }
        }
    }
    std::cout<<"Total points: "<<pts.size()<<std::endl;

    auto test_map_ptr = std::make_shared<TestMap>(filename);
    PathPreprocessor<TestMap> preprocessor(test_map_ptr,pts);
    preprocessor.process();

}

int main(int argc, char **argv) {


    //path_preprocessor_test();
    bgi_distance_test();
    //planner_test();

    /*
    std::string half_track_data_file = "../data/tracks/temp_nwh_half.csv";
    auto reader = CSVReader(half_track_data_file);
    auto raw_data = reader.read();
    DM waypoints;
    if(!raw_data.toDM(waypoints)){
        std::cout<<"Read Track File Fails\n";
    }
    std::cout<<waypoints<<std::endl;
    acsr::Path path(waypoints,7,50);
    path.plot();*/

    //track_test();




    //::testing::InitGoogleTest(&argc, argv);
    //return RUN_ALL_TESTS();
    return 0;

}
