//
// Created by lab on 8/31/22.
//

#ifndef RACECARPLANNER_TIRE_MODEL_HPP
#define RACECARPLANNER_TIRE_MODEL_HPP
#include <casadi/casadi.hpp>

namespace acsr{

    class PacejkaModel {
        using param_t = std::map<std::string,double>;
    public:
        PacejkaModel() = default;
        PacejkaModel(const param_t& param){
            B_long = param.at("B_long");
            C_long = param.at("C_long");
            D_long = param.at("D_long");
            B_lat = param.at("B_lat");
            C_lat = param.at("C_lat");
            D_lat = param.at("D_lat");
        };

        void get_force(){

        }

        double get_lambda_at_max_force(){
            return tan(M_PI/(2*C_long))/B_long;
        }

        double get_alpha_at_max_force(){
            return tan(M_PI/(2*C_lat))/B_lat;
        }

        template<class T>
        T get_max_longitudinal_force(T& Fz){
            return Fz*D_long;
        }



    private:
        double B_long,C_long,D_long,B_lat,C_lat,D_lat;
    };



    class PacejkaSimpleModel {
        using param_t = std::map<std::string,double>;
    public:
        PacejkaSimpleModel() = default;
        PacejkaSimpleModel(const param_t& param){
            B_long = param.at("B_long");
            C_long = param.at("C_long");
            D_long = param.at("D_long");
            B_lat = param.at("B_lat");
            C_lat = param.at("C_lat");
            D_lat = param.at("D_lat");
        }

        template<class T>
        T get_force(T& lamb,T& alpha,double Fz){
            return T::vertcat(Fz*D_long*casadi::sin(C_long*casadi::atan(B_long*lamb)),
                                 Fz*D_lat*casadi::sin(C_lat*casadi::atan(B_lat*alpha)));
        }

        template<class T>
        T get_longitudinal_force(T& lamb,double Fz) {
            return Fz * D_long * casadi::sin(C_long * casadi::atan(B_long * lamb));
        }

        template<class T>
        T get_lateral_force(T& alpha,double Fz) {
            return Fz * D_lat * casadi::sin(C_lat * casadi::atan(B_lat * alpha));
        }


        double get_lambda_at_max_force(){
            return tan(M_PI/(2*C_long))/B_long;
        }

        double get_alpha_at_max_force(){
            return tan(M_PI/(2*C_lat))/B_lat;
        }

        template<class T>
        T get_max_longitudinal_force(T& Fz){
            return Fz*D_long;
        }

    private:
        double B_long,C_long,D_long,B_lat,C_lat,D_lat;
    };
}

#endif //RACECARPLANNER_TIRE_MODEL_HPP
