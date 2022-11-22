//
// Created by acsr on 11/17/22.

#include "global_planner_optimizer.hpp"
#include <matplotlibcpp.h>
#include "nlohmann/json.hpp"

using namespace acsr;

int main(){

    std::string database_file{"../output/global_planner.db"};
    DM waypoints;
    CSVReader reader("../data/temp.csv");
    auto csv_data = reader.read();
    if(!csv_data.toDM(waypoints)){
        std::cout<<"Convert csv data to dm fail\n";
        return 1;
    }

    std::cout<<waypoints.rows()<<std::endl;
    auto path_ptr = std::make_shared<Path>(waypoints);


    std::string vehicle_file = "../data/params/racecar.json";
    if(!std::filesystem::exists(vehicle_file)){
        std::cout<<vehicle_file<<" does not exist\n";
        return 1;
    }
    std::ifstream ifs(vehicle_file);
    auto vehicle_params = json::parse(ifs);
    double track_width = 10.0;

    int N = 3000;

    BicycleKinematicOptimizer optimizer(path_ptr,track_width,vehicle_params);
    param_t init{{"s",0.0},{"v",0.15}};

    auto return_value = optimizer.make_plan(init,-1,N);



    if(return_value.first){
        DM dt = std::get<0>(return_value.second);
        DM x = std::get<1>(return_value.second);
        DM u = std::get<2>(return_value.second);


        SQLite::Database    db(database_file, SQLite::OPEN_READWRITE|SQLite::OPEN_CREATE);
        optimizer.save(db, dt,x,u,std::vector<std::string>{"t","n","phi","v"}, std::vector<std::string>{"delta","d"});

        std::cout<<x(Slice(0,3),Slice())<<std::endl;

        namespace plt = matplotlibcpp;
        std::vector<double> dt_vec(dt->size()+1);
        std::partial_sum(dt->begin(),dt->end(),dt_vec.begin()+1);
        dt_vec[0]=0.0;

        auto tau = DM::linspace(0,path_ptr->get_max_tau(),500+1).T();
        auto n = DM::zeros(1,500+1);
        auto n_ones = DM::ones(1,500+1);
        auto centerline = path_ptr->f_tn_to_xy(DMVector{tau,n})[0];
        auto dm_centerline_x = centerline(0,Slice());
        auto dm_centerline_y = centerline(1,Slice());

        auto innerline = path_ptr->f_tn_to_xy(DMVector{tau,n_ones*track_width/2})[0];
        auto dm_innerline_x = innerline(0,Slice());
        auto dm_innerline_y = innerline(1,Slice());

        auto outerline = path_ptr->f_tn_to_xy(DMVector{tau,-n_ones*track_width/2})[0];
        auto dm_outerline_x = outerline(0,Slice());
        auto dm_outerline_y = outerline(1,Slice());

        plt::plot(std::vector<double>(dm_centerline_x->begin(),dm_centerline_x->end()),std::vector<double>(dm_centerline_y->begin(),dm_centerline_y->end()),"g--",
                  std::vector<double>(dm_innerline_x->begin(),dm_innerline_x->end()),std::vector<double>(dm_innerline_y->begin(),dm_innerline_y->end()),"b-",
                  std::vector<double>(dm_outerline_x->begin(),dm_outerline_x->end()),std::vector<double>(dm_outerline_y->begin(),dm_outerline_y->end()),"b-");

        auto path_line = path_ptr->f_tn_to_xy(DMVector{x(0,Slice()),x(1,Slice())})[0];
        auto path_line_x = path_line(0,Slice());
        auto path_line_y = path_line(1,Slice());
        plt::plot(std::vector<double>(path_line_x->begin(),path_line_x->end()),std::vector<double>(path_line_y->begin(),path_line_y->end()),"r-");

        plt::show();

    }

}

//
