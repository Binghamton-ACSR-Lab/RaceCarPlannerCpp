//
// Created by lab on 9/1/22.
//

#ifndef RACECARPLANNER_AXLE_HPP
#define RACECARPLANNER_AXLE_HPP

#include "wheel.hpp"
namespace acsr {
    template<typename ScalarType, typename VectorType>
    class axle_t {
        axle_t()=default;
    private:
        wheel<ScalarType,VectorType> left_wheel;
        wheel<ScalarType,VectorType> right_wheel;
    };
}

#endif //RACECARPLANNER_AXLE_HPP
