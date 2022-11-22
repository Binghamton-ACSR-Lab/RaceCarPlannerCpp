//
// Created by acsr on 11/17/22.



#include "websocketpp/server.hpp"
#include "websocketpp/config/asio_no_tls.hpp"
#include "nlohmann/json.hpp"
#include "websocket_server.hpp"
#include <fstream>


int main(){
    acsr::WsServer server;
    server.run([](const json& value){
        std::vector<double> x = value["waypoints"]["ptsx"];
        std::vector<double> y = value["waypoints"]["ptsy"];
        std::ofstream outfile("../data/temp.csv");
        for(auto i=0;i<x.size();i++)
            outfile<<std::fixed << std::setprecision(1)<<std::setw(8)<<x[i]<<", "<<std::setw(10)<<y[i]<<std::endl;
        outfile.close();
        },
               [](const json& value){
        std::cout<<value<<std::endl;

        return R"({"foo": "bar"})"_json;
    });
}


//
