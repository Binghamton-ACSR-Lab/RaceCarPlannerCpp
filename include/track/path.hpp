//
// Created by acsr on 9/19/22.
//

#ifndef RACECARPLANNER_PATH_HPP
#define RACECARPLANNER_PATH_HPP


#include <string>
#include <iostream>
#include "rtree.hpp"
#include "casadi/casadi.hpp"
//#include <matplotlibcpp.h>


using casadi::Slice;
using namespace casadi;

namespace acsr {
    /***
     * a class representing a spline path based on CasAdi
     */
    class Path {
    public:
        Path() = delete;

        explicit Path(DM pts,bool closed=false) {
            waypoints_=pts;
            if(waypoints_.rows()!=2){
                assert(pts.columns()==2);
                waypoints_ = pts.T();
            }else{
                waypoints_=pts;
            }

            f_xy = get_parametric_function(waypoints_);
            t_max = waypoints_.columns()-(closed?0:1);


            auto ts = casadi::DM::linspace(0, t_max, resolution * t_max);
            auto ts_vec = std::vector<double>(ts->begin(),ts->end());

            auto t = MX::sym("t");
            auto n = MX::sym("n");
            auto s = MX::sym("s");
            auto pt_t_mx = f_xy(t)[0];
            auto jac = MX::jacobian(pt_t_mx,t);
            auto hes=MX::jacobian(jac,t);

            f_tangent_vec = casadi::Function("vec",{t},{jac});
            auto kappa = (jac(0)*hes(1)-jac(1)*hes(0))/MX::pow(MX::norm_2(jac),3);
            f_kappa = casadi::Function("kappa",{t},{kappa});
            f_phi = casadi::Function("phi",{t},{casadi::MX::atan2(jac(1),jac(0))});

            const auto theta = M_PI/2;
            auto rot_mat_vector = std::vector<std::vector<double>>{{cos(theta),-sin(theta)},{sin(theta),cos(theta)}};
            auto rot_mat = DM(rot_mat_vector);
            MXDict dae;
            dae["x"] = s;
            dae["t"] = t;
            dae["ode"] = MX::norm_2(f_tangent_vec(t)[0]);

            Function integ = integrator("inte","cvodes",dae,{{"grid",ts_vec}});
            DMDict args;
            args["x0"] =  0;
            auto s_inte = integ(args);
            auto s_value = s_inte["xf"];
            auto s_value_vec = std::vector<double>(s_value->begin(),s_value->end());
            s_value_vec.insert(s_value_vec.begin(),0);
            s_max = s_value_vec.back();

            s_to_t_lookup = interpolant("s_to_t","linear",std::vector<std::vector<double>>{s_value_vec},ts_vec);
            t_to_s_lookup = interpolant("t_to_s","linear",std::vector<std::vector<double>>{ts_vec},s_value_vec);

            auto xy =pt_t_mx + MX::mtimes(rot_mat,jac/MX::norm_2(jac))*n;

            f_tn_to_xy = casadi::Function("tn_to_xy",{t,n},{xy});

            auto phi = f_phi(t)[0];
            auto dm_xy = MX::sym("dm_xy",2,1);
            auto dm_n = -MX::sin(phi)*(dm_xy-pt_t_mx)(0) + MX::cos(phi)*(dm_xy-pt_t_mx)(1);
            f_xy_to_tn = casadi::Function("xy_to_tn",{dm_xy,t},{dm_n});

            DM ns_zeros = DM::zeros(ts.rows(),ts.columns());
            auto center_line = f_tn_to_xy(DMVector {ts.T(),ns_zeros.T()})[0];
            auto dm_x = center_line(0,Slice());
            auto dm_y = center_line(1,Slice());
            std::vector<double> vec_x{dm_x->begin(),dm_x->end()};
            std::vector<double> vec_y{dm_y->begin(),dm_y->end()};
            std::vector<double> vec_ts{ts->begin(),ts->end()};

            search_tree = RTree(vec_x,vec_y,vec_ts);

        }

        DM xy2t(const DM& xy,bool refined=false){
            DM dm_x,dm_y;

            if(xy.is_vector()){
                dm_x=xy(0);
                dm_y=xy(1);
            }else {
                dm_x = xy(0, Slice());
                dm_y = xy(1, Slice());
            }
            std::vector<double> vec_x{dm_x->begin(),dm_x->end()};
            std::vector<double> vec_y{dm_y->begin(),dm_y->end()};
            auto ts_vec = search_tree.findNearest(vec_x,vec_y);
            //return DM(ts_vec);
            if(!refined){
                return DM(ts_vec);
            }
            //auto tangent_vec = f_tangent_vec(dm_t)[0];
            //auto xy_on_line = f_xy(dm_t)[0];

            DM tau_real(1,vec_x.size());
            for(auto i=0;i<vec_x.size();++i) {
                auto t = MX::sym("t");
                auto tangent_vec = f_tangent_vec(t)[0];
                auto xy_on_line = f_xy(t)[0];

                auto norm_vec = xy(Slice(),i) - xy_on_line;
                auto dot_prod = norm_vec(0) * tangent_vec(0) + norm_vec(1) * tangent_vec(1);
                auto f = casadi::Function("f", {t}, {dot_prod});
                auto rf = casadi::rootfinder("rf", "newton", f);
                auto t_real = rf(DM(ts_vec[i]))[0];
                tau_real(0,i)=t_real(0);
            }

            //DM n = f_xy_to_tn(DMVector{xy,t_real})[0];
            return tau_real;

        }




/*
        void plot(double width = 10.0,int resolution=100){
            namespace plt = matplotlibcpp;
            auto ts = casadi::DM::linspace(0, t_max, resolution * t_max);
            DM ns_zeros = DM::zeros(ts.rows(),ts.columns());
            DM ns_ones = DM::ones(ts.rows(),ts.columns());
            auto center_line = f_tn_to_xy(DMVector {ts.T(),ns_zeros.T()})[0];
            auto inner_line = f_tn_to_xy(std::vector<DM>{ts.T(),+width/2*ns_ones.T()})[0];
            auto outer_line = f_tn_to_xy(std::vector<DM>{ts.T(),-width/2*ns_ones.T()})[0];
            auto center_x_dm = center_line(0,Slice());
            auto center_y_dm = center_line(1,Slice());
            auto inner_x_dm = inner_line(0,Slice());
            auto inner_y_dm = inner_line(1,Slice());
            auto outer_x_dm = outer_line(0,Slice());
            auto outer_y_dm = outer_line(1,Slice());

            plt::plot(std::vector<double>{center_x_dm->begin(),center_x_dm->end()},
                      std::vector<double>{center_y_dm->begin(),center_y_dm->end()},
                      "r-",
                      std::vector<double>{inner_x_dm->begin(),inner_x_dm->end()},
                      std::vector<double>{inner_y_dm->begin(),inner_y_dm->end()},
                      "k-",
                      std::vector<double>{outer_x_dm->begin(),outer_x_dm->end()},
                      std::vector<double>{outer_y_dm->begin(),outer_y_dm->end()},
                      "b-");

            auto pts_x_dm = waypoints_(0,Slice());
            auto pts_y_dm = waypoints_(1,Slice());

            plt::scatter(std::vector<double>{pts_x_dm->begin(),pts_x_dm->end()},
                      std::vector<double>{pts_y_dm->begin(),pts_y_dm->end()},10);
            plt::show();
        }
*/
        double get_max_length(){
            return s_max;
        }

        double get_max_tau(){
            return t_max;
        }

    public:
        /**
         * casadi function, return (x,y) from tau, tau is the paramter of the parametric spline
         */
        casadi::Function f_xy;

        /**
         * casadi function, return curvature from tau
         */
        casadi::Function f_kappa;

        /**
         * casadi function, return heading of the spline at current tau
         */
        casadi::Function f_phi;

        /**
         * return a vector(dx/dt,dy/dt) of the spline at t
         */
        casadi::Function f_tangent_vec;

        /**
         * convert length to tau
         */
        casadi::Function s_to_t_lookup;

        /**
         * convert tau to length
         */
        casadi::Function t_to_s_lookup;

        /**
         * convert point in (t,n) frame to (x,y) frame
         */
        casadi::Function f_tn_to_xy;

        /**
         * get n in (t,n) frame of a point in (x,y) frame. two parameters of this function, (xy,tau)
         */
        casadi::Function f_xy_to_tn;
    private:
        casadi::DM waypoints_;
        double t_max;
        double s_max;

        const int resolution = 100;

        //casadi::DM center_line;
        //casadi::DM inner_line;
        //casadi::DM outer_line;
        RTree search_tree;


        static casadi::Function get_parametric_function(const casadi::DM& waypoints,bool closed= false){

            assert(waypoints.rows()==2);
            auto t = casadi::MX::sym("t");
            auto dm_x = waypoints(0,Slice());
            auto dm_y = waypoints(1,Slice());
            auto size = waypoints.columns()+(closed?1:0);
            std::vector<double> tau_vec(size),x(size),y(size);
            std::generate(tau_vec.begin(), tau_vec.end(), [i = 0]() mutable {return i++;});
            std::copy(dm_x->begin(),dm_x->end(),x.begin());
            std::copy(dm_y->begin(),dm_y->end(),y.begin());
            if(closed){
                x.back() = x.front();
                y.back() = y.front();
            }

            auto fx =casadi::interpolant("f_x","bspline",std::vector<std::vector<double>>{tau_vec},x);
            auto fy = casadi::interpolant("f_y","bspline",std::vector<std::vector<double>>{tau_vec},y);

            auto g = casadi::MX::vertcat({fx(t)[0],fy(t)[0]});
            casadi::Function f_pt_t("f", casadi::MXVector{t}, casadi::MXVector{g}, {"t"}, {"pt"});
            return f_pt_t;

        }
    };
}


#endif //RACECARPLANNER_PATH_HPP
