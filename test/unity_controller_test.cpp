//
// Created by acsr on 11/21/22.
//

#include "unity_controller.hpp"

using namespace acsr;
int main(){
    UnityController<acsr::BicycleDynamicsTwoBrakeOptimizer<PacejkaSimpleModel>> controller;
    controller.run();

}
