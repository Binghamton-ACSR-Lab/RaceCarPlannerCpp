//
// Created by acsr on 9/22/22.
//

#ifndef RACECARPLANNER_TEST_MAP_HPP
#define RACECARPLANNER_TEST_MAP_HPP

#include <tinyxml2.h>
#include <string>
#include <iostream>
#include "shape.hpp"
#include <boost/geometry.hpp>

namespace bg = boost::geometry;
typedef bg::model::point<double, 2, bg::cs::cartesian> point_t;
typedef bg::model::polygon<point_t> polygon_t;

namespace acsr{
    class TestMap{
    public:
        TestMap()=default;

        explicit TestMap(const std::string& filename,double ratio){
            tinyxml2::XMLDocument doc;
            doc.LoadFile( filename.c_str());
            auto root = doc.RootElement();
            for (tinyxml2::XMLElement* child = root->FirstChildElement(); child != nullptr; child = child->NextSiblingElement())
            {
                if(strcmp(child->Name(),"rectangle")==0){
                    auto width = child->DoubleAttribute("width")*ratio;
                    auto height = child->DoubleAttribute("height")*ratio;
                    auto x = child->DoubleAttribute("x")*ratio;
                    auto y = child->DoubleAttribute("y")*ratio;
                    add_rectangle(x,y,width,height);
                }else if(strcmp(child->Name(),"circle")==0){
                    auto x = child->DoubleAttribute("x")*ratio;
                    auto y = child->DoubleAttribute("y")*ratio;
                    auto radius = child->DoubleAttribute("radius")*ratio;
                    add_circle(x,y,radius);
                }

            }
            std::cout<<shapes.size()<<'\n';
        }

        bool valid(double x,double y){
            return std::all_of(std::begin(shapes),std::end(shapes),[x,y](auto &s){
                return !bg::within(point_t {x,y},s);
            });
        }

        template<class T>
        bool valid(T& state){
            return valid(state[0],state[1]);
        }


        template<class T>
        std::vector<std::pair<double,T>> near_collides(double distance_holder,double x,double y){
            point_t pt{x,y};
            std::vector<std::pair<double,T>> result;
            //DistancePoint<point_t> distance_point{};
            //distance_point.distance = std::numeric_limits<double>::max();
            for(auto& s : shapes) {
                auto dist_poly = point_to_polygon(pt, s);

                if (dist_poly.distance < distance_holder) {
                    result.push_back(std::make_pair(dist_poly.distance,T{dist_poly.projected_point.get<0>(),dist_poly.projected_point.get<1>()}));
                }
            }
            return result;
        }

        template<class T>
        std::vector<T> get_force(std::vector<T> const& vec,double max_edge_margin){

            std::vector<T> result;
            for(auto i=0;i<vec.size();++i){
                point_t pt{vec[i][0],vec[i][1]};
                T force{0,0};

                for(auto& s : shapes) {
                    auto dist_poly = point_to_polygon(pt, s);

                    if (dist_poly.distance < max_edge_margin) {
                        auto dx = (pt.get<0>() - dist_poly.projected_point.get<0>()) / dist_poly.distance;
                        auto dy = (pt.get<1>() - dist_poly.projected_point.get<1>()) / dist_poly.distance;
                        auto k = force_ratio/(dist_poly.distance*dist_poly.distance);
                        if(bg::within(point_t {pt.get<0>(),pt.get<1>()},s))
                            k=-k;
                        force = T{force[0] + k * dx, force[1] +k*dy};
                    }
                }
                result.push_back(force);

            }
            return result;
        }

        template<class T>
        std::vector<T> get_force(std::vector<T> const& vec){

            std::vector<T> result;
            for(auto i=0;i<vec.size();++i){
                point_t pt{vec[i][0],vec[i][1]};
                T force{0,0};

                for(auto& s : shapes) {
                    auto dist_poly = point_to_polygon(pt, s);

                    //if (dist_poly.distance < max_edge_margin) {
                    auto dx = (pt.get<0>() - dist_poly.projected_point.get<0>()) / dist_poly.distance;
                    auto dy = (pt.get<1>() - dist_poly.projected_point.get<1>()) / dist_poly.distance;
                    auto k = force_ratio/(dist_poly.distance*dist_poly.distance);
                    if(bg::within(point_t {pt.get<0>(),pt.get<1>()},s))
                        k=-k;
                    force = T{force[0] + k * dx, force[1] +k*dy};
                    //}
                }
                result.push_back(force);

            }
            return result;
        }

        std::vector<std::pair<std::vector<double>,std::vector<double>>> plot_data(){
            std::vector<std::pair<std::vector<double>,std::vector<double>>> pts;
            for(auto& shape: shapes) {
                std::vector<double> x,y;
                bg::for_each_point(shape,[&x,&y](point_t& pt){
                    x.push_back(pt.get<0>());
                    y.push_back(pt.get<1>());
                });
                pts.emplace_back(x,y);
            }
            return pts;
            /*
            std::vector<std::pair<std::vector<double>,std::vector<double>>> v(shapes.size());
            std::transform(shapes.begin(),shapes.end(),v.begin(),[](point_t& s){
                return s->get_plot_data();
            });
            return v;*/
        }

        void add_circle(double x,double y,double radius){
            polygon_t poly;
            for(auto i=0;i<101;++i){
                bg::append(poly.outer(), point_t(x+radius*cos(2*M_PI/100*i), y+radius*sin(2*M_PI/100*i)));
            }
            shapes.push_back(poly);
        }

        void add_rectangle(double x,double y,double width,double height){
            polygon_t poly;
            bg::append(poly.outer(), point_t(x, y));
            bg::append(poly.outer(), point_t(x+width, y));
            bg::append(poly.outer(), point_t(x+width, y+height));
            bg::append(poly.outer(), point_t(x, y+height));
            bg::append(poly.outer(), point_t(x, y));
            shapes.push_back(poly);
        }

        void add_polygon(std::initializer_list<point_t> pts){
            shapes.push_back(polygon_t{pts});
        }



    private:
        //std::vector<std::shared_ptr<Shape<double>>> shapes;
        //std::vector<Circle<double>> circles;
        std::vector<polygon_t> shapes;
        const double force_ratio = 1e3;
    };

}

#endif //RACECARPLANNER_TEST_MAP_HPP
