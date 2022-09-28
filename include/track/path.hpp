//
// Created by acsr on 9/19/22.
//

#ifndef RACECARPLANNER_PATH_HPP
#define RACECARPLANNER_PATH_HPP


#include <string>
#include <iostream>
#include "rtree.hpp"
#include <matplotlibcpp.h>


using casadi::Slice;
using namespace casadi;
namespace acsr {

    class Path {
    public:
        Path() = delete;

        explicit Path(DM pts, double max_width,int resolution = 100):max_width_(max_width) {
            waypoints_=pts;
            if(waypoints_.rows()!=2){
                assert(pts.columns()==2);
                waypoints_ = pts.T();
            }else{
                waypoints_=pts;
            }

            pt_t = get_parametric_function(waypoints_);
            t_max = waypoints_.columns()-1;

            auto ts = casadi::DM::linspace(0, t_max, resolution * t_max);
            auto ts_vec = std::vector<double>(ts->begin(),ts->end());

            auto t = MX::sym("t");
            auto n = MX::sym("n");
            auto s = MX::sym("s");
            auto pt_t_mx = pt_t(t)[0];
            auto jac = MX::jacobian(pt_t_mx,t);
            auto hes=MX::jacobian(jac,t);

            f_tangent_vec = casadi::Function("vec",{t},{jac});
            auto kappa = (jac(0)*hes(1)-jac(1)*hes(0))/MX::pow(MX::norm_2(jac),3);
            f_kappa = casadi::Function("kappa",{t},{kappa});
            f_phi = casadi::Function("phi",{t},{casadi::MX::atan2(jac(1),jac(0))});

            //auto test_phi = f_phi(ts.T())[0];
            //std::cout<<test_phi.size()<<std::endl;
            //std::cout<<test_phi<<std::endl;

            auto theta = M_PI/2;
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


            DM ns_zeros = DM::zeros(ts.rows(),ts.columns());
            DM ns_ones = DM::ones(ts.rows(),ts.columns());

            auto xy =pt_t_mx + MX::mtimes(rot_mat,jac/MX::norm_2(jac))*n;

            f_tn_to_xy = casadi::Function("tn_to_xy",{t,n},{xy});

            auto phi = f_phi(t)[0];
            //auto phi_rot_mat_vector = MX(2,2);
            //phi_rot_mat_vector(0,0) = MX::cos(phi);
            //phi_rot_mat_vector(0,1) = MX::sin(phi);
            //phi_rot_mat_vector(1,0) = -MX::sin(phi);
            //phi_rot_mat_vector(1,1) = MX::cos(phi);

            MX dm_xy;
            auto dm_n = -MX::sin(phi)*(dm_xy-pt_t_mx)(0) + MX::cos(phi)*(dm_xy-pt_t_mx)(1);
            //auto tn = MX::mtimes(phi_rot_mat_vector,dm_xy-pt_t_mx);
            f_xy_to_tn = casadi::Function("xy_to_tn",{dm_xy,t},{dm_n});
            //MX::mtimes(phi_rot_mat_vector,)

            center_line = f_tn_to_xy(DMVector {ts.T(),ns_zeros.T()})[0];
            //std::cout<<center_line_test.size()<<std::endl;

            inner_line = f_tn_to_xy(std::vector<DM>{ts.T(),+max_width_/2*ns_ones.T()})[0];
            outer_line = f_tn_to_xy(std::vector<DM>{ts.T(),-max_width_/2*ns_ones.T()})[0];

            //std::cout<<casadi::DM::sum2(casadi::DM::sum1(center_line.T()-center_line_test))<<std::endl;
            //std::cout<<casadi::DM::sum2(casadi::DM::sum1(outer_line.T()-outer_line_test))<<std::endl;
            //std::cout<<casadi::DM::sum2(casadi::DM::sum1(inner_line.T()-inner_line_test))<<std::endl;

            search_tree = RTree(center_line);

        }



        void plot(){
            namespace plt = matplotlibcpp;
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

        double get_max_length(){
            return s_max;
        }

        double get_max_tau(){
            return t_max;
        }

        constexpr double get_width(){
            return max_width_;
        }

    public:
        casadi::Function pt_t;
        casadi::Function f_kappa;
        casadi::Function f_phi;
        casadi::Function f_tangent_vec;
        casadi::Function s_to_t_lookup;
        casadi::Function t_to_s_lookup;
        casadi::Function f_tn_to_xy;
        casadi::Function f_xy_to_tn;
    private:
        casadi::DM waypoints_;
        const double max_width_;

        double t_max;
        double s_max;
        casadi::DM center_line;
        casadi::DM inner_line;
        casadi::DM outer_line;

        RTree search_tree;


        static casadi::Function get_parametric_function(const casadi::DM& waypoints){

            assert(waypoints.rows()==2);
            auto t = casadi::MX::sym("t");

            std::vector<double> tau_vec(waypoints.columns());
            auto dm_x = waypoints(0,Slice());
            auto dm_y = waypoints(1,Slice());
            auto x = std::vector<double>{dm_x->begin(),dm_x->end()};
            auto y = std::vector<double>{dm_y->begin(),dm_y->end()};
            std::generate(tau_vec.begin(), tau_vec.end(), [i = 0]() mutable {return i++;});

            auto fx =casadi::interpolant("f_x","bspline",std::vector<std::vector<double>>{tau_vec},x);
            auto fy = casadi::interpolant("f_y","bspline",std::vector<std::vector<double>>{tau_vec},y);

            auto g = casadi::MX::vertcat({fx(t)[0],fy(t)[0]});
            casadi::Function f_pt_t("f", casadi::MXVector{t}, casadi::MXVector{g}, {"t"}, {"pt"});
            return f_pt_t;

        }
    };
}


#endif //RACECARPLANNER_PATH_HPP
