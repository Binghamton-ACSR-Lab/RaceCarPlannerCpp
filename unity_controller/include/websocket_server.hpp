//
// Created by acsr on 11/8/22.
//

#ifndef CAR_CONTROL_CPP_WEBSOCKET_SERVER_H
#define CAR_CONTROL_CPP_WEBSOCKET_SERVER_H

#include "websocketpp/server.hpp"
#include "websocketpp/config/asio_no_tls.hpp"
#include "nlohmann/json.hpp"
#include <fstream>
#include "chrono"

using json = nlohmann::json;
using namespace std::chrono_literals;

using ws_server_t = websocketpp::server<websocketpp::config::asio>;
namespace acsr{
    class WsServer{
    public:
        WsServer(double dt):dt_(dt) {}


        void run(std::function<json(json)> f1,std::function<json(json)> f2){
            server.set_message_handler(std::bind(&WsServer::on_message<std::function<json(json)>,std::function<json(json)>>, this,std::placeholders::_1, std::placeholders::_2,f1,f2));
            server.set_access_channels(websocketpp::log::alevel::all);
            server.set_error_channels(websocketpp::log::elevel::all);

            server.init_asio();
            server.listen(4567);
            server.start_accept();
            server.run();
        }

    private:
        ws_server_t server;
        std::chrono::duration<double> dt_;

        template<class F1,class F2>
        void on_message(websocketpp::connection_hdl hdl, ws_server_t::message_ptr msg,F1 f1,F2 f2){
            auto start = std::chrono::system_clock::now();

            auto msg_str = msg->get_payload();
            if(msg_str.substr(0,2) != "42"){
                std::cout<<"message type error\n";
                return;
            }
            auto json_obj = json::parse(msg_str.substr(2));
            if(json_obj.empty() || json_obj[0]!="telemetry")return;

            json& data = json_obj[1];
            if(data.contains("waypoints")){
                json s = f1(data);
                server.send(hdl, R"(42["reference",)" + s.dump() + "]", msg->get_opcode());
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
                // Some computation here
                auto end = std::chrono::system_clock::now();
                std::chrono::duration<double> elapsed_seconds = end-start;
                if(elapsed_seconds<dt_-10ms){
                    std::this_thread::sleep_for(dt_-10ms-elapsed_seconds);
                }
                end = std::chrono::system_clock::now();
                std::chrono::duration<double> elapsed = end-start;
                std::cout << "Waited " << elapsed.count() << " s\n";

                server.send(hdl, R"(42["steer",)" + s.dump() + "]", msg->get_opcode());
            }
        }
    };



}

#endif //CAR_CONTROL_CPP_WEBSOCKET_SERVER_H
