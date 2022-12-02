//
// Created by xli185 on 9/25/22.
//

#ifndef RACECARPLANNER_TEST_H
#define RACECARPLANNER_TEST_H

#include <boost/geometry.hpp>
#include <iostream>
#include "utility"
#include "test_map.hpp"
#include "waypoints_processor.hpp"
#include "global_planner_optimizer.hpp"

using namespace acsr;

void bgi_distance_test(){
    namespace bg = boost::geometry;
    typedef bg::model::point<double, 2, bg::cs::cartesian> point_t;
    typedef bg::model::polygon<point_t> polygon_t;
    polygon_t polygon{{{0.0, 0.0}, {0.0, 5.0}, {5.0, 5.0}, {5.0, 0.0}, {0.0, 0.0}}};

    point_t pt{2,8};


    auto segments = boost::make_iterator_range(bg::segments_begin(polygon), bg::segments_end(polygon));
    for(auto segment : segments){
        //std::cout<<"("<<segment.first->get<0>()<<","<<segment.first->get<1>()<<")\t"
        //                                                                  "("<<segment.second->get<0>()<<","<<segment.second->get<1>()<<")\n";

        auto pt_with_distance = point_to_segment(pt,*segment.first,*segment.second);
        std::cout<<pt_with_distance.distance<<": "<<pt_with_distance.projected_point.get<0>()<<","<<pt_with_distance.projected_point.get<1>()<<'\n';
    }

    auto com_distance = bg::comparable_distance(polygon,pt);
    std::cout<<double(com_distance);


}

/*
void make_planner_test(){
    namespace plt = matplotlibcpp;

    //read the waypoints
    std::vector<std::vector<double>> pts;
    auto map_file = std::string("../data/map/test_map.xml");
    tinyxml2::XMLDocument doc;
    doc.LoadFile( map_file.c_str());
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

    //create the map, this map is used as valid_checker.
    auto test_map_ptr = std::make_shared<TestMap>(map_file);

    double edge_margin = 13.0;

    // create the waypoint processor, used to smoothing waypoints
    //template valid_checker_t has to be implemented two functions:
    // bool valid(x,y) & std::vector<T> get_force(std::vector<T> const& vec)
    WaypointsProcessor<TestMap> processor(test_map_ptr,pts);

    //get the smoothed waypoints
    auto refined_waypoints = processor.process();
    std::cout<<refined_waypoints<<std::endl;

    //create a spline path going through those smoothed waypoints
    std::shared_ptr<Path> path_ptr = std::make_shared<Path>(refined_waypoints);

    //load params
    std::string vehicle_file = "../data/params/racecar_nwh.yaml";
    std::string front_tire_file = "../data/params/nwh_tire.yaml";
    std::string rear_tire_file = "../data/params/nwh_tire.yaml";

    if(!std::filesystem::exists(vehicle_file)){
        std::cout<<vehicle_file<<" does not exist\n";
        return;
    }
    if(!std::filesystem::exists(front_tire_file)){
        std::cout<<front_tire_file<<" does not exist\n";
        return;
    }
    if(!std::filesystem::exists(rear_tire_file)){
        std::cout<<front_tire_file<<" does not exist\n";
        return;
    }
    param_t vehicle_params,front_tire_params,rear_tire_params;
    auto vehicle_yaml = YAML::LoadFile(vehicle_file);
    std::cout<<"load vehicle config file.. Total node: "<<vehicle_yaml.size()<<std::endl;
    for(auto it = vehicle_yaml.begin();it!=vehicle_yaml.end();++it){
        vehicle_params[it->first.as<std::string>()]=it->second.as<double>();
    }
    auto front_tire_yaml = YAML::LoadFile(front_tire_file);
    std::cout<<"load front tire config file.. Total node: "<<front_tire_yaml.size()<<std::endl;
    for(auto it = front_tire_yaml.begin();it!=front_tire_yaml.end();++it){
        front_tire_params[it->first.as<std::string>()]=it->second.as<double>();
    }
    if(front_tire_file==rear_tire_file){
        rear_tire_params = front_tire_params;
    }else {
        auto rear_tire_yaml = YAML::LoadFile(rear_tire_file);
        std::cout<<"load rear tire config file.. Total node: "<<rear_tire_yaml.size()<<std::endl;
        for (auto it = rear_tire_yaml.begin(); it != rear_tire_yaml.end(); ++it) {
            rear_tire_params[it->first.as<std::string>()] = it->second.as<double>();
        }
    }

    //create the global planner optimizer
    //BicycleDynamicsTwoBrakeOptimizer<TestMap> optimizer(path_ptr,2*edge_margin,vehicle_params,front_tire_params,rear_tire_params,150,test_map_ptr);
    //get the trajectory
    //auto optimized_path = optimizer.make_plan();

    //plot
    plt::figure(1);
    auto collision = test_map_ptr->plot_data();
    for(auto& data:collision){
        plt::plot(data.first,data.second,"k-");
    }

    DM original_waypoints = DM(pts).T();
    auto original_x = original_waypoints(0,Slice());
    auto original_y = original_waypoints(1,Slice());
    plt::scatter(std::vector<double>(original_x->begin(),original_x->end()),std::vector<double>(original_y->begin(),original_y->end()),30.0,{{"c", "red"}, {"marker","*"}});

    auto refined_x = refined_waypoints(0,Slice());
    auto refined_y = refined_waypoints(1,Slice());
    plt::scatter(std::vector<double>(refined_x->begin(),refined_x->end()),std::vector<double>(refined_y->begin(),refined_y->end()),30.0,{{"c", "blue"}, {"marker","o"}});

    auto s = DM::linspace(0,path_ptr->get_max_length(),201).T();
    auto t = path_ptr->s_to_t_lookup(s)[0];
    auto n = DM::zeros(1,201);

    auto path_xy = path_ptr->f_tn_to_xy(DMVector{t,n})[0];
    auto path_x = path_xy(0,Slice());
    auto path_y = path_xy(1,Slice());
    plt::named_plot("path",std::vector<double>(path_x->begin(),path_x->end()),std::vector<double>(path_y->begin(),path_y->end()));

    auto boundary_xy = optimizer.boundary_xy();
    auto outer_x = boundary_xy.first(0,Slice());
    auto outer_y = boundary_xy.first(1,Slice());
    auto inner_x = boundary_xy.second(0,Slice());
    auto inner_y = boundary_xy.second(1,Slice());
    plt::named_plot("outer boundary",std::vector<double>(outer_x->begin(),outer_x->end()),std::vector<double>(outer_y->begin(),outer_y->end()));
    plt::named_plot("inner boundary",std::vector<double>(inner_x->begin(),inner_x->end()),std::vector<double>(inner_y->begin(),inner_y->end()));


    if(optimized_path.first){
        auto x = std::get<1>(optimized_path.second);
        auto dm_xy = path_ptr->f_tn_to_xy(std::vector<DM>{x(0,Slice()),x(1,Slice())})[0];
        auto trajectory_x = dm_xy(0,Slice());
        auto trajectory_y = dm_xy(1,Slice());
        plt::named_plot("trajectory",std::vector<double>(trajectory_x->begin(),trajectory_x->end()),std::vector<double>(trajectory_y->begin(),trajectory_y->end()));
        plt::legend();

        plt::figure(2);
        auto dt = std::get<0>(optimized_path.second);
        auto vx = x(3,Slice(0,-1));
        auto vy = x(4,Slice(0,-1));

        std::vector<double> t(dt.columns());
        t[0]=0;
        std::partial_sum(dt->begin(),dt->end()-1,t.begin()+1);

        plt::named_plot("vx",std::vector<double>(t.begin(),t.end()),std::vector<double>(vx->begin(),vx->end()));
        plt::named_plot("vy",std::vector<double>(t.begin(),t.end()),std::vector<double>(vy->begin(),vy->end()));

    }

    plt::legend();
    plt::show();
}
*/

#endif //RACECARPLANNER_TEST_H
