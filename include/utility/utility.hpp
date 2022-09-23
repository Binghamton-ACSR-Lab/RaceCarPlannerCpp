//
// Created by lab on 8/31/22.
//

#ifndef RACECARPLANNER_UTILITY_HPP
#define RACECARPLANNER_UTILITY_HPP

#include <casadi/casadi.hpp>

namespace acsr{
    template<typename T>
    void printArray(T& arr){
        for(auto& t:arr){
            std::cout<<t<<'\t';
        }
        std::cout<<'\n';
    }
}

#endif //RACECARPLANNER_UTILITY_HPP
