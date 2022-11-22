//
// Created by acsr on 11/8/22.
//

#ifndef CAR_CONTROL_CPP_WEBSOCKET_SERVER_H
#define CAR_CONTROL_CPP_WEBSOCKET_SERVER_H

#include "websocketpp/server.hpp"
#include "websocketpp/config/asio_no_tls.hpp"
#include "nlohmann/json.hpp"
#include <fstream>

using json = nlohmann::json;


using ws_server_t = websocketpp::server<websocketpp::config::asio>;
namespace acsr{
    class WsServer{
    public:
        WsServer() = default;

        void run(std::function<void(json)> f1,std::function<json(json)> f2){
            server.set_message_handler(std::bind(&WsServer::on_message<std::function<void(json)>,std::function<json(json)>>, this,std::placeholders::_1, std::placeholders::_2,f1,f2));
            server.set_access_channels(websocketpp::log::alevel::all);
            server.set_error_channels(websocketpp::log::elevel::all);

            server.init_asio();
            server.listen(4567);
            server.start_accept();
            server.run();
        }

    private:
        ws_server_t server;

        template<class F1,class F2>
        void on_message(websocketpp::connection_hdl hdl, ws_server_t::message_ptr msg,F1 f1,F2 f2){
            auto msg_str = msg->get_payload();
            if(msg_str.substr(0,2) != "42"){
                std::cout<<"message type error\n";
                return;
            }
            auto json_obj = json::parse(msg_str.substr(2));
            if(json_obj.empty() || json_obj[0]!="telemetry")return;

            json& data = json_obj[1];
            if(data.contains("waypoints")){
                f1(data);
                /*
                std::vector<double> x = data["waypoints"]["ptsx"];
                std::vector<double> y = data["waypoints"]["ptsy"];
                std::ofstream outfile("../data/temp.csv");
                for(auto i=0;i<x.size();i++)
                    outfile<<std::fixed << std::setprecision(1)<<std::setw(8)<<x[i]<<", "<<std::setw(10)<<y[i]<<std::endl;
                outfile.close();
                return;*/
            }else {
                json s = f2(data);
                server.send(hdl, R"(42["steer",)" + s.dump() + "]", msg->get_opcode());
            }
        }
    };



}

#endif //CAR_CONTROL_CPP_WEBSOCKET_SERVER_H
