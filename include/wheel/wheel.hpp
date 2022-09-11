//
// Created by lab on 8/31/22.
//

#ifndef RACECARPLANNER_WHEEL_HPP
#define RACECARPLANNER_WHEEL_HPP

#include "tire_model.hpp"

namespace acsr{

    template<typename ScalarType, typename VectorType>
    class wheel{

    public:
        wheel(){

        }

        //! Calls Tire::update(x0,v0,omega) of the base class, and calls update_self()
        //! This updates the tire dynamics and tire forces
        //! @param[in] x0: new frame origin position [m]
        //! @param[in] v0: new frame origin velocity [m/s]
        //! @param[in] omega: new value for tire angular speed [rad/s]
        void update(cosnt VectorType& x0,  const VectorType& v0, ScalarType omega){

        }

    private:
        tire_model_t tire_model_;

    };



}


#endif //RACECARPLANNER_WHEEL_HPP
