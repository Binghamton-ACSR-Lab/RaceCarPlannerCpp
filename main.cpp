#include <iostream>

#include <casadi/casadi.hpp>
#include "include/utility/utility.hpp"
#include <gtest/gtest.h>
#include "track.hpp"
using namespace acsr;


template<class T>
T update(T& x,T& u){
    return T::vertcat({x,u});
}

int main(int argc, char **argv) {
    //std::vector<std::vector<double>> vec1{{1,2,3,4},{5,6,7,8}};
    casadi::DM vv1{1,3,4,5};
    casadi::DM vv2{2,3,4,5};
    std::cout<<casadi::DM::vertcat({vv1.T(),vv2.T()})<<std::endl;
    //auto vv2 = v;

    Track track("../data/tracks/temp.csv",5,true);
    track.plot();
    //SX tf =SX::sym("tf");
    std::vector<double> tf{10.0};
    //double tf = 10.0;
    int nx = 3;
    int np = 1;
    SX x  = SX::sym("x",nx);  // state
    SX p  = SX::sym("u",np);  // control

    // ODE right hand side function
    SX ode = vertcat((1 - x(1)*x(1))*x(0) - x(1) + p,
                     x(0),
                     x(0)*x(0) + x(1)*x(1) + p*p);
    SXDict dae = {{"x", x}, {"p", p}, {"ode", ode}};
    Function ref_integrator = integrator("ref_integrator",
                                         "cvodes", dae, {{"grid",tf},{"output_t0",true}});

    DMDict args;
    args["x0"] =   DM{1.0,2.0,3.0};
    auto s_inte = ref_integrator(args);
    std::cout<<s_inte;


    std::vector<std::vector<double>> vec{{1,2},{3,4},{5,6},{7,8},{9,10}};
    DM dm(vec);

    auto x_axis = dm(Slice(),0);
    auto y_axis = dm(Slice(),1);

    std::vector<double> x_vec = std::vector<double>{x_axis->begin(),x_axis->end()};
    std::vector<double> y_vec = std::vector<double>{y_axis->begin(),y_axis->end()};

    for(auto i=0;i<x_vec.size();++i)
        std::cout<<x_vec[i]<<'\t';
    std::cout<<std::endl;

    casadi::MX t(1,1);
    auto s = 2*t+2;

    DM a{0};


    casadi::DM b{1};
    auto c = update(a,b);
    std::cout<<c<<std::endl;

    casadi::MX aa(1);
    casadi::MX bb(1);

    auto cc = update(aa,bb);
    std::cout<<cc<<std::endl;

    auto f = casadi::Function("f",{t},{s});



    std::cout<<f(DM{2})<<std::endl;


    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();


}
