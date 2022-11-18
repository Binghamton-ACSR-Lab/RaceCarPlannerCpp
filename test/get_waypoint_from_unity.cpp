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
        std::cout<<value<<std::endl;

        return R"({"foo": "bar"})"_json;
    });
}


//
