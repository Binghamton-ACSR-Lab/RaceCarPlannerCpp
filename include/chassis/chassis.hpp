//
// Created by lab on 9/1/22.
//

#ifndef RACECARPLANNER_CHASSIS_HPP
#define RACECARPLANNER_CHASSIS_HPP

#include "axle.hpp"
namespace acsr{
    template<typename ScalarType, typename VectorType>
    class chassis{
    public:
        chassis() = default;
    private:
        axle_t<ScalarType,VectorType> front_axle_;
        axle_t<ScalarType,VectorType> rear_axle_;

        ScalarType u_; //! [in] Longitudinal velocity (in road frame) [m/s]
        ScalarType v_; //! [in] Lateral velocity (in road frame) [m/s]
        ScalarType omega_; //! [in] Yaw speed [rad/s]

        double m_; //mass
        double Iz_; //inertia

        double x_front_axle_;
        double x_rear_axle_;


    };
}

#endif //RACECARPLANNER_CHASSIS_HPP
