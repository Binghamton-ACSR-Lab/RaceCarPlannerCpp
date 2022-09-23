//
// Created by acsr on 9/22/22.
//

#ifndef RACECARPLANNER_TEST_MAP_HPP
#define RACECARPLANNER_TEST_MAP_HPP

#include <tinyxml2.h>
#include <string>
#include <iostream>
#include "shape.hpp"

namespace acsr{
    class TestMap{
    public:
        TestMap()=default;
        TestMap(const std::string& filename){
            tinyxml2::XMLDocument doc;
            doc.LoadFile( filename.c_str());
            auto root = doc.RootElement();
            for (tinyxml2::XMLElement* child = root->FirstChildElement(); child != nullptr; child = child->NextSiblingElement())
            {
                if(strcmp(child->Name(),"rectangle")==0){
                    auto width = child->DoubleAttribute("width");
                    auto height = child->DoubleAttribute("height");
                    auto x = child->DoubleAttribute("x");
                    auto y = child->DoubleAttribute("y");
                    shapes.push_back(std::make_shared<Rectangle<double>>(x,y,width,height));
                }else if(strcmp(child->Name(),"circle")==0){
                    auto x = child->DoubleAttribute("x");
                    auto y = child->DoubleAttribute("y");
                    auto radius = child->DoubleAttribute("radius");
                    shapes.push_back(std::make_shared<Circle<double>>(x,y,radius));
                }

            }
            std::cout<<shapes.size()<<'\n';
            for(auto& s:shapes){
                std::cout<<s->type()<<", Area: "<<s->area()<<'\n';
            }
        }

        bool valid(double x,double y){
            return std::all_of(std::begin(shapes),std::end(shapes),[x,y](auto &s){
                return s->collision(x,y)==false;
            });
        }

        std::vector<std::pair<std::vector<double>,std::vector<double>>> plot_data(){
            std::vector<std::pair<std::vector<double>,std::vector<double>>> v(shapes.size());
            std::transform(shapes.begin(),shapes.end(),v.begin(),[](auto& s){
                return s->get_plot_data();
            });
            return v;
        }

    private:
        std::vector<std::shared_ptr<Shape<double>>> shapes;
        //std::vector<Circle<double>> circles;

    };

}

#endif //RACECARPLANNER_TEST_MAP_HPP
