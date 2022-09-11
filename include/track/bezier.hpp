//
// Created by lab on 9/2/22.
//

#ifndef RACECARPLANNER_BEZIER_HPP
#define RACECARPLANNER_BEZIER_HPP
#include <string>
#include <iostream>
#include "csv_reader.hpp"

namespace acsr {


    class CubicBezier {
    public:
        //CubicBezier() = default;
        static std::tuple<casadi::DM,casadi::DM> get_coef(const casadi::DM& waypoints){
            using namespace casadi;
            auto all = Slice();
            auto n = waypoints.rows()-1;
            auto M = casadi::DM::zeros(n,n);
            const auto tridiagle = casadi::DM{1, 4, 1};
            for(auto idx =0;idx<n-2;++idx){
                M(idx+1,casadi::Slice(idx,idx+3)) = tridiagle;
            }
            M(0,Slice(0,2)) = tridiagle(Slice(1,3));
            M(0,n-1) = 1;
            M(n-1,Slice(n-2,n))= tridiagle(Slice(0,2));
            M(n-1,0)= 1;

            auto s =DM::zeros(n,2);
            for(auto idx =0;idx<n;++idx){
                s(idx,all) = 2*(2*waypoints(idx,all) + waypoints(idx+1,all));
            }

            auto Ax = DM::solve(M,s(all,0));
            auto Ay = DM::solve(M,s(all,1));
            auto a = DM::vertcat({Ax,Ay});
            a = DM::reshape(a,n,2);
            auto b = DM::zeros(n,2);
            b(Slice(0,n-1),all) = 2* waypoints(Slice(1,n),all) - a(Slice(1,n),all);
            b(n-1,all) = 2*waypoints(0,all) - a(0,all);
            return std::tie(a,b);

        }




    };
}

#endif //RACECARPLANNER_TRACK_HPP
