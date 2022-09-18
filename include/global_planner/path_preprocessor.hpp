//
// Created by acsr on 9/18/22.
//

#ifndef RACECARPLANNER_PATH_PREPROCESSOR_HPP
#define RACECARPLANNER_PATH_PREPROCESSOR_HPP

#include "casadi/casadi.hpp"
#include <future>

using namespace casadi;
using param_t = std::map<std::string,double>;
template<class map_t>
class PathPreprocessor{

public:
    PathPreprocessor(map_t* map, DM& pts,const param_t& vehicle_params, int segments = 2,double min_distance = 1.0,double min_edge_margin = 1.0);
    PathPreprocessor(map_t* map, const std::vector<std::vector<double>>& pts,const param_t& vehicle_params,int segments = 2,double min_distance = 1.0,double min_edge_margin = 1.0);

    DM operator()(){

    }

private:
    int segments_;
    double min_distance_;
    double min_edge_margin_;

    std::vector<int> split(){
        if(segments_==1)return std::vector<int>{};

        std::vector<int> result(segments_-1);

    }


};


#endif //RACECARPLANNER_PATH_PREPROCESSOR_HPP
