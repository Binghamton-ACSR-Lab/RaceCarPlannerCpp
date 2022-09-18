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
        static std::tuple<casadi::DM,casadi::DM> get_coef(const casadi::DM& waypoints,bool closed){
            using namespace casadi;
            auto all = Slice();
            const auto tridiagle = casadi::DM{1, 4, 1};
            casadi_int n;
            DM M,s;
            if(closed){
                n = waypoints.rows();
                M = casadi::DM::zeros(n,n);
                s =DM::zeros(n,2);
                for(auto idx =0;idx<n-2;++idx){
                    M(idx+1,casadi::Slice(idx,idx+3)) = tridiagle;
                }
                M(0,Slice(0,2)) = tridiagle(Slice(1,3));
                M(0,n-1) = 1;
                M(n-1,Slice(n-2,n))= tridiagle(Slice(0,2));
                M(n-1,0)= 1;
                for(auto idx =0;idx<n-1;++idx){
                    s(idx,all) = 2*(2*waypoints(idx,all) + waypoints(idx+1,all));
                }
                s(n-1,all) = 2*(2*waypoints(n-1,all) + waypoints(0,all));
            }else{
                n = waypoints.rows()-1;
                M = casadi::DM::zeros(n,n);
                s =DM::zeros(n,2);
                for(auto idx =1;idx<n-2;++idx){
                    M(idx+1,casadi::Slice(idx,idx+3)) = tridiagle;
                }
                M(0,Slice(0,2)) = DM{2,1};
                M(n-1,Slice(n-2,n))= DM{2,7};
                for(auto idx =1;idx<n-1;++idx){
                    s(idx,all) = 2*(2*waypoints(idx,all) + waypoints(idx+1,all));
                }
                s(0,all) = 2*waypoints(0,all) + waypoints(1,all);
                s(n-1,all) = 8*waypoints(n-1,all) + waypoints(n,all);
            }

            auto Ax = DM::solve(M,s(all,0));
            auto Ay = DM::solve(M,s(all,1));
            auto a = DM::vertcat({Ax,Ay});
            a = DM::reshape(a,n,2);
            auto b = DM::zeros(n,2);
            b(Slice(0,n-1),all) = 2* waypoints(Slice(1,n),all) - a(Slice(1,n),all);
            if(closed)
                b(n-1,all) = 2*waypoints(0,all) - a(0,all);
            else
                b(n-1,all) = (waypoints(n,all) + a(n-1,all))/2;
            return std::tie(a,b);

        }


        /*
         *
        # find the a & b points
        def get_bezier_coef(points):
            # since the formulas work given that we have n+1 points
            # then n must be this:
            n = len(points) - 1

            # build coefficents matrix
            C = 4 * np.identity(n)
            np.fill_diagonal(C[1:], 1)
            np.fill_diagonal(C[:, 1:], 1)
            C[0, 0] = 2
            C[n - 1, n - 1] = 7
            C[n - 1, n - 2] = 2

            # build points vector
            P = [2 * (2 * points[i] + points[i + 1]) for i in range(n)]
            P[0] = points[0] + 2 * points[1]
            P[n - 1] = 8 * points[n - 1] + points[n]

            # solve system, find a & b
            A = np.linalg.solve(C, P)
            B = [0] * n
            for i in range(n - 1):
                B[i] = 2 * points[i + 1] - A[i + 1]
            B[n - 1] = (A[n - 1] + points[n]) / 2

            return A, B

        # returns the general Bezier cubic formula given 4 control points
        def get_cubic(a, b, c, d):
            return lambda t: np.power(1 - t, 3) * a + 3 * np.power(1 - t, 2) * t * b + 3 * (1 - t) * np.power(t, 2) * c + np.power(t, 3) * d

         */


    };
}

#endif //RACECARPLANNER_TRACK_HPP
