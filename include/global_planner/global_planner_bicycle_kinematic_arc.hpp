//
// Created by acsr on 11/16/22.
//

#ifndef RACECARPLANNER_GLOBAL_PLANNER_BICYCLE_KINEMATIC_ARC_HPP
#define RACECARPLANNER_GLOBAL_PLANNER_BICYCLE_KINEMATIC_ARC_HPP

#include "utility.hpp"
#include "path.hpp"

namespace acsr {

    template<class valid_checker_t=void>
    class BicycleKinematicOptimizer {
    public:
        BicycleKinematicOptimizer() = default;

        explicit BicycleKinematicOptimizer(std::shared_ptr<Path> path_ptr,
                                           double path_width,
                                           const param_t &vehicle_params,
                                           int optimization_resolution = 100,
                                           std::shared_ptr<valid_checker_t> valid_checker = nullptr)
                : vehicle_params_(vehicle_params), path_width_(path_width), N_(optimization_resolution),
                  path_(path_ptr) {

            option_["max_iter"] = 3000;
            option_["tol"] = 1e-6;
            //option_["linear_solver"]="ma57";

            if (std::is_void<valid_checker_t>::value || valid_checker == nullptr) {
                road_constraint_upper_bound_ = DM{path_width / 2};
                road_constraint_lower_bound_ = DM{-path_width / 2};
            } else {
                std::vector<size_t> index_vect(N_ + 1);
                std::iota(index_vect.begin(), index_vect.end(), 0);

                auto s = DM::linspace(0, path_->get_max_length(), N_ + 1);
                s = s.T();
                auto t = path_->s_to_t_lookup(s)[0];
                auto n = DM::zeros(1, N_ + 1);
                auto center_line = path_->f_tn_to_xy(DMVector{t, n})[0];

                road_constraint_upper_bound_ = path_width / 2 * DM::ones(1, N_ + 1);
                road_constraint_lower_bound_ = -path_width / 2 * DM::ones(1, N_ + 1);

                std::for_each(index_vect.begin(), index_vect.end(), [&](size_t idx) {
                    std::vector<std::pair<double, DM>> result = valid_checker->template near_collides<DM>(
                            path_width / 2, double(center_line(0, idx)), double(center_line(1, idx)));
                    for (auto &r: result) {
                        auto n = path_->f_xy_to_tn(DMVector{r.second, t(0, idx)})[0](0);
                        if (double(n(0)) < 0 && -r.first > double(road_constraint_lower_bound_(0, idx))) {
                            road_constraint_lower_bound_(0, idx) = -r.first;
                        } else if (double(n(0)) > 0 && r.first < double(road_constraint_upper_bound_(0, idx))) {
                            road_constraint_upper_bound_(0, idx) = r.first;
                        }
                    }
                });
                std::cout << road_constraint_lower_bound_ << std::endl;
                std::cout << road_constraint_upper_bound_ << std::endl;


                if (path_width / 2 - double(DM::mmin(road_constraint_upper_bound_)) < 0.1 &&
                    path_width / 2 + double(DM::mmax(road_constraint_lower_bound_)) < 0.1) {
                    road_constraint_upper_bound_ = DM{path_width / 2};
                    road_constraint_lower_bound_ = DM{-path_width / 2};
                }

            }

        }

        std::pair<bool, std::tuple<DM, DM, DM>> make_plan(double n0 = 0, double v0 = 0.15, bool print = true) {


            auto all = Slice();
            auto _0_N = Slice(0, N_);
            auto _N_1 = Slice(1, N_ + 1);

            Dict casadi_option;
            casadi_option["print_time"] = print;
            if (print) {
                option_["print_level"] = 5;
            } else {
                option_["print_level"] = 1;
            }

            //const auto track_width = path_->get_width();
            const auto delta_min = vehicle_params_.at("delta_min");
            const auto delta_max = vehicle_params_.at("delta_max");
            const auto v_min = vehicle_params_.at("v_min");
            const auto v_max = vehicle_params_.at("v_max");
            const auto a_min = vehicle_params_.at("a_min");
            const auto a_max = vehicle_params_.at("a_max");


            const auto s_val = casadi::DM::linspace(0, path_->get_max_length(), N_ + 1).T();
            const auto tau_array = path_->s_to_t_lookup(s_val)[0];
            const auto tau0 = double(tau_array(0));
            const auto phi0 = double(path_->f_phi(tau_array(0))[0]);
            const auto tau_array_mid = path_->s_to_t_lookup((s_val(0, _0_N) + s_val(0, _N_1)) * 0.5);

            const auto kappa_array = path_->f_kappa(tau_array_mid)[0];
            //const auto kappa_array = (kappa_array_temp(0,_N_1) + kappa_array_temp(0,Slice(1,N+1)))*0.5;
            const auto tangent_vec_array = path_->f_tangent_vec(tau_array_mid)[0];
            //const auto tangent_vec_array = (tangent_vec_array_temp(all,_N_1) + tangent_vec_array_temp(all,Slice(1,N+1)))*0.5;

            const auto phi_array = DM::atan2(tangent_vec_array(1, all), tangent_vec_array(0, all));
            const auto tangent_vec_norm = DM::sqrt((tangent_vec_array(0, all) * tangent_vec_array(0, all) +
                                                    tangent_vec_array(1, all) * tangent_vec_array(1, all)));

            //tau_array.to_file("/home/acsr/Documents/data/tau_array_cpp.txt");
            //kappa_array.to_file("/home/acsr/Documents/data/kappa_array_cpp.txt");
            //tangent_vec_array.to_file("/home/acsr/Documents/data/tangent_vec_array_cpp.txt");
            //phi_array.to_file("/home/acsr/Documents/data/phi_array_cpp.txt");
            //tangent_vec_norm.to_file("/home/acsr/Documents/data/tangent_vec_norm_cpp.txt");

            casadi::Opti opti;
            auto X = opti.variable(nx, N_ + 1);
            auto U = opti.variable(nu, N_);
            //auto X_dot = opti.variable(nx, N);
            auto dt_sym_array = opti.variable(1, N_);

            auto n_sym_array = (X(IDX_X_n, _0_N) + X(IDX_X_n, _N_1)) * 0.5;
            auto phi_sym_array = (X(IDX_X_phi, _0_N) + X(IDX_X_phi, _N_1)) * 0.5;
            auto vx_sym_array = (X(IDX_X_vx, _0_N) + X(IDX_X_vx, _N_1)) * 0.5;
            //auto omega_sym_array = (X(IDX_X_omega,_0_N)+X(IDX_X_omega,_N_1))*0.5;

            auto delta_sym_array = U(IDX_U_delta, all);
            auto a_sym_array = U(IDX_U_a, all);

            //auto n_sym_array = X(1,Slice());
            //auto n_obj = (casadi::MX::atan(5 * ( n_sym_array*n_sym_array - track->get_width()*track->get_width() / 4) ) + casadi::pi / 2) * 12;

            auto X0 = casadi::DM::vertcat({tau0, n0, phi0, v0});
            auto dphi_c_sym_array = phi_sym_array - phi_array;

            //X_dot(0, all) = dt_sym_array*(vx_sym_array * MX::cos(dphi_c_sym_array) - vy_sym_array * MX::sin(dphi_c_sym_array))/(tangent_vec_norm*(1-n_sym_array*kappa_array));
            //X_dot(1, all) = dt_sym_array*(vx_sym_array * MX::sin(dphi_c_sym_array) + vy_sym_array* MX::cos(dphi_c_sym_array)) ;
            //X_dot(2, all) = dt_sym_array*omega_sym_array;
            //X_dot(3, all) = dt_sym_array * (fx_r_sym_array + fx_f_sym_array * MX::cos(delta_sym_array) - fy_f_sym_array * MX::sin(delta_sym_array) + mass * vy_sym_array* omega_sym_array)/ mass;
            //X_dot(4, all) = dt_sym_array * (fy_r_sym_array + fx_f_sym_array * MX::sin(delta_sym_array) + fy_f_sym_array * MX::cos(delta_sym_array) - mass * vx_sym_array * omega_sym_array)/ mass;
            //X_dot(5, all) = dt_sym_array * (fy_f_sym_array * lf * MX::cos(delta_sym_array) + fx_f_sym_array * lf * MX::sin(delta_sym_array) - fy_r_sym_array * lr)/ Iz ;
            //X_dot(6, all) = dt_sym_array * delta_dot_sym_array;

            auto upbound = DM::mmin(road_constraint_upper_bound_);
            auto lowbound = DM::mmax(road_constraint_lower_bound_);

            //auto n_low_obj = MX::exp(MX(road_constraint_lower_bound_)-X(IDX_X_n,all));
            //auto n_upper_obj = MX::exp(X(IDX_X_n,all)-MX(road_constraint_upper_bound_));
            // auto n_low_obj = MX::atan(10*(MX(road_constraint_lower_bound_)-X(IDX_X_n,all)));
            // auto n_upper_obj = MX::atan(10*(X(IDX_X_n,all)-MX(road_constraint_upper_bound_)));

            //auto n_obj = MX::atan(5 * (n_sym_array * n_sym_array - (path_width_ / 2) *(path_width_ / 2))) + casadi::pi / 2;
            //opti.minimize(MX::sum2(dt_sym_array) + MX::dot(n_obj, n_obj));
            //opti.minimize(MX::sum2(dt_sym_array) + MX::dot(delta_dot_sym_array,delta_dot_sym_array) + 15.0*MX::dot(n_obj,n_obj));
            //opti.minimize(MX::sum2(dt_sym_array) + MX::dot(delta_dot_sym_array,delta_dot_sym_array));
            opti.minimize(10 * MX::sum2(
                    dt_sym_array)); /*+ MX::dot(delta_dot_sym_array,delta_dot_sym_array) + MX::sum2(n_low_obj) + MX::sum2(n_upper_obj));*/
            //dynamics
            //opti.subject_to(X(all, Slice(1,N+1)) == X(all, _N_1) + X_dot);
            opti.subject_to(X(0, _N_1) == X(0, _0_N) + dt_sym_array * (vx_sym_array * MX::cos(dphi_c_sym_array)) /
                                                       (tangent_vec_norm * (1 - n_sym_array * kappa_array)));
            opti.subject_to(X(1, _N_1) == X(1, _0_N) + dt_sym_array * (vx_sym_array * MX::sin(dphi_c_sym_array)));
            //opti.subject_to(X(2, _N_1) == X(2, _0_N) + dt_sym_array * vx_sym_array * MX::sin(dphi_c_sym_array)/vehicle_params_["wheel_base"]);
            opti.subject_to(X(2, _N_1) == X(2, _0_N) + dt_sym_array * vx_sym_array * MX::tan(delta_sym_array) /
                                                       vehicle_params_.at("wheel_base"));
            opti.subject_to(X(3, _N_1) == X(3, _0_N) + dt_sym_array * a_sym_array);

            // opti.subject_to(opti.bounded(-0.1,X(1, -1),0.1));
            //inital conditions
            opti.subject_to(X(0, all) == tau_array);


            opti.subject_to(X(all, 0) == X0);
            opti.subject_to(dt_sym_array > 0);


            //state boundary
            opti.subject_to(opti.bounded(delta_min, delta_sym_array, delta_max));
            //opti.subject_to(opti.bounded(road_constraint_lower_bound_,X(IDX_X_n,all),road_constraint_upper_bound_));
            //opti.subject_to(opti.bounded(delta_min, delta_sym_array, delta_max));
            opti.subject_to(opti.bounded(lowbound, n_sym_array, upbound));
            //control boundary
            //opti.subject_to(opti.bounded(delta_dot_min, delta_dot_sym_array, delta_dot_max));
            opti.subject_to(opti.bounded(a_min, a_sym_array, a_max));
            opti.subject_to(opti.bounded(v_min, vx_sym_array, v_max));
            // opti.subject_to(opti.bounded(road_constraint_lower_bound_,n_sym_array,road_constraint_upper_bound_));

            auto X_guess = casadi::DM::zeros(nx, N_ + 1);
            //X_guess(0,Slice()) = tau_array;
            //X_guess(2,Slice(1,N+1)) = phi_array;
            //X_guess(3,all) = v0;
            for (auto i = 0; i < N_ + 1; ++i) {
                X_guess(Slice(), i) = X0;
            }
            opti.set_initial(X, X_guess);
            opti.set_initial(dt_sym_array, 0.1);


            opti.solver("ipopt", casadi_option, option_);
            try {
                auto sol = opti.solve();
                auto dt_array = sol.value(dt_sym_array);
                auto sol_x = sol.value(X);
                auto sol_u = sol.value(U);
                return std::make_pair(true, std::make_tuple(dt_array, sol_x, sol_u));
            }
            catch (CasadiException &e) {
                std::cout << e.what() << std::endl;
                std::cout << "Solve Optimal Problem Fails\n";
                return std::make_pair(false, std::make_tuple(DM{}, DM{}, DM{}));;
            }

        }

        std::pair<DM, DM> boundary_xy() {
            auto s = DM::linspace(0, path_->get_max_length(), N_ + 1);
            s = s.T();
            auto t = path_->s_to_t_lookup(s)[0];
            //auto n = DM::zeros(1,N_+1);
            auto outer_line = path_->f_tn_to_xy(DMVector{t, road_constraint_lower_bound_})[0];
            auto inner_line = path_->f_tn_to_xy(DMVector{t, road_constraint_upper_bound_})[0];
            return std::make_pair(outer_line, inner_line);
        }

        std::pair<DM, DM> boundary_n() {
            return std::make_pair(road_constraint_lower_bound_, road_constraint_upper_bound_);
        }


    private:
        std::shared_ptr<Path> path_;
        double path_width_{};
        param_t vehicle_params_, track_params_;

        Dict option_;
        DM road_constraint_upper_bound_, road_constraint_lower_bound_;
        const int N_{};

        //std::shared_ptr<valid_checker_t> valid_checker_;

        constexpr static int nx = 4;
        constexpr static int nu = 2;

        constexpr static unsigned IDX_X_t = 0;
        constexpr static unsigned IDX_X_n = 1;
        constexpr static unsigned IDX_X_phi = 2;
        constexpr static unsigned IDX_X_vx = 3;

        constexpr static unsigned IDX_U_delta = 0;
        constexpr static unsigned IDX_U_a = 1;


    };

}


#endif //RACECARPLANNER_GLOBAL_PLANNER_BICYCLE_KINEMATIC_ARC_HPP
