//
// Created by acsr on 9/2/22.
//

#ifndef RACECARPLANNER_TRACK_HPP
#define RACECARPLANNER_TRACK_HPP
#include <string>
#include <iostream>
#include "csv_reader.hpp"
#include "bezier.hpp"
#include "rtree.hpp"
#include <matplotlibcpp.h>

using casadi::Slice;
using namespace casadi;
namespace acsr {

    class Track {
    public:
        Track() = delete;

        Track(const std::string &file_name, double track_width,bool is_closed,int resolution = 100) {
            auto reader = CSVReader(file_name);
            auto raw_data = reader.read();
            if(!raw_data.toDM(waypoints)){
                std::cout<<"Read Track File Fails\n";
                return;
            }
            auto rows = waypoints.rows();
            auto cols = waypoints.columns();
            if(is_closed){
                waypoints.resize(rows+1,cols);
                waypoints(rows,Slice()) = waypoints(0,Slice());
                t_max = rows;
            }else{
                t_max = rows-1;
            }
            auto ts = casadi::DM::linspace(0, t_max, resolution * t_max);
            auto ts_vec = std::vector<double>(ts->begin(),ts->end());



            //dm_a.to_file("a.txt");
            //dm_b.to_file("b.txt");

            pt_t = get_parametric_function(waypoints);

            auto t = MX::sym("t");
            auto s = MX::sym("s");
            auto pt_t_mx = pt_t(t)[0];
            auto jac = MX::jacobian(pt_t_mx,t);
            auto hes=MX::jacobian(jac,t);

            f_tangent_vec = casadi::Function("vec",{t},{jac});

            auto kappa = (jac(0)*hes(1)-jac(1)*hes(0))/MX::pow(MX::norm_2(jac),3);
            f_kappa = casadi::Function("kappa",{t},{kappa});

            //auto vec_mx = f_tangent_vec(t)[0];
            //f_tangent_vec = casadi::Function("phi",{t},{vec});
            f_phi = casadi::Function("phi",{t},{casadi::MX::atan2(jac(1),jac(0))});

            auto test_phi = f_phi(ts.T())[0];
            std::cout<<test_phi.size()<<std::endl;
            std::cout<<test_phi<<std::endl;

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

            //center_line.resize(resolution*t_max,2);
            //inner_line.resize(resolution*t_max,2);
            //outer_line.resize(resolution*t_max,2);

            center_line = pt_t(ts.T())[0];
            //std::cout<<ts.size()<<std::endl;
            //std::cout<<center_line.size()<<std::endl;
            center_line = casadi::DM::reshape(center_line,2,7800);

            //std::cout<<center_line(Slice(),0)<<std::endl;
            //std::cout<<center_line(Slice(),1)<<std::endl;
            //std::cout<<waypoints(0,Slice())<<std::endl;
            //std::cout<<waypoints(1,Slice())<<std::endl;


            auto vec = f_tangent_vec(ts.T())[0];
            //std::cout<<vec.size()<<std::endl;
            //vec = casadi::DM::reshape(vec,2,7800);
            auto norm = casadi::DM::sqrt(vec(0,Slice())*vec(0,Slice())+vec(1,Slice())*vec(1,Slice()));
            //std::cout<<norm.size()<<std::endl;
            vec = vec/casadi::DM::repmat(norm,2);
            //std::cout<<vec.size()<<std::endl;
            vec = DM::mtimes(rot_mat,vec);
            inner_line = (center_line+track_width/2*vec).T();
            outer_line = (center_line-track_width/2*vec).T();
            center_line = center_line.T();
            /*
            for(auto idx=0;idx<resolution*t_max;++idx){
                auto v = pt_t(ts(idx))[0];
                v = DM::reshape(v,1,2);

                center_line(idx,Slice()) = v;
                auto vec = f_tangent_vec(ts)[0];
                vec = vec/casadi::DM::norm_2(vec);
                vec = DM::mtimes(rot_mat,vec);
                inner_line(idx,Slice()) = v+track_width/2*vec.T();
                outer_line(idx,Slice()) = v-track_width/2*vec.T();
            }*/
            search_tree = RTree(center_line);

        }



        void plot(){
            namespace plt = matplotlibcpp;
            auto center_x_dm = center_line(Slice(),0);
            auto center_y_dm = center_line(Slice(),1);
            auto inner_x_dm = inner_line(Slice(),0);
            auto inner_y_dm = inner_line(Slice(),1);
            auto outer_x_dm = outer_line(Slice(),0);
            auto outer_y_dm = outer_line(Slice(),1);

            plt::plot(std::vector<double>{center_x_dm->begin(),center_x_dm->end()},
                      std::vector<double>{center_y_dm->begin(),center_y_dm->end()},
                      "r-",
                      std::vector<double>{inner_x_dm->begin(),inner_x_dm->end()},
                      std::vector<double>{inner_y_dm->begin(),inner_y_dm->end()},
                      "k-",
                      std::vector<double>{outer_x_dm->begin(),outer_x_dm->end()},
                      std::vector<double>{outer_y_dm->begin(),outer_y_dm->end()},
                      "b-");

            std::vector<double> inner_dot_x,inner_dot_y,outer_dot_x,outer_dot_y,center_dot_x,center_dot_y;
            for(auto i=0;i<center_line.rows();i+=100){
                inner_dot_x.push_back(double(inner_line(i,0)));
                inner_dot_y.push_back(double(inner_line(i,1)));

                outer_dot_x.push_back(double(outer_line(i,0)));
                outer_dot_y.push_back(double(outer_line(i,1)));

                center_dot_x.push_back(double(center_line(i,0)));
                center_dot_y.push_back(double(center_line(i,1)));
            }
            plt::scatter(inner_dot_x,inner_dot_y);
            plt::scatter(outer_dot_x,outer_dot_y);
            plt::scatter(center_dot_x,center_dot_y);


            plt::show();
        }

        double getMaxLength(){
            return s_max;
        }

    public:
        casadi::Function pt_t;
        casadi::Function f_kappa;
        casadi::Function f_phi;
        casadi::Function f_tangent_vec;
        casadi::Function s_to_t_lookup;
        casadi::Function t_to_s_lookup;

    private:
        casadi::DM waypoints;

        double t_max;
        double s_max;
        casadi::DM center_line;
        casadi::DM inner_line;
        casadi::DM outer_line;

        RTree search_tree;


        casadi::Function get_parametric_function(const casadi::DM& waypoints){
            auto [dm_a,dm_b] = CubicBezier::get_coef(waypoints);
            auto dm_A = casadi::MX(dm_a);
            auto dm_B = casadi::MX(dm_b);
            auto dm_waypoints = casadi::MX(waypoints);
            auto t = casadi::MX::sym("t");
            auto n = dm_a.rows();
            assert(dm_b.rows()==n);
            auto tau = casadi::MX::mod(t,n);

            auto i = casadi::MX::floor(tau);
            auto a = dm_waypoints(i,Slice());

            auto b = dm_A(i,Slice());

            auto c = dm_B(i,Slice());
            auto i1 = casadi::MX::mod(i+1,n);
            auto d =dm_waypoints(i1,Slice());
            auto g = casadi::MX::pow(1 - (tau-i), 3) * a + 3 * casadi::MX::pow(1 - (tau-i), 2) * (tau-i) * b
                    + 3 * (1 - (tau-i)) * casadi::MX::pow(tau-i, 2) * c + casadi::MX::pow(tau-i, 3) * d;
            //auto vec = casadi::MX::jacobian(g,tau);

            casadi::Function f_pt_t("f",casadi::MXVector{t},casadi::MXVector{g},{"t"},{"pt"});
            //casadi::Function f_tangent_vec("vec",casadi::MXVector{t},casadi::MXVector{vec},{"t"},{"vec"});
            //return std::tuple(f_pt_t,f_tangent_vec);
            return f_pt_t;
        }
    };
}

#endif //RACECARPLANNER_TRACK_HPP
