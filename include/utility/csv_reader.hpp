//
// Created by lab on 9/2/22.
//

#ifndef RACECARPLANNER_CSV_READER_HPP
#define RACECARPLANNER_CSV_READER_HPP

#include <string>
#include <vector>
#include <fstream>
#include <boost/algorithm/string.hpp>
#include <casadi/casadi.hpp>
#include <charconv>
#include <boost/lexical_cast.hpp>
#include <regex>

namespace acsr {
    class CSVData {
    public:
        CSVData(std::vector<std::string> header,std::vector<std::vector<std::string>> data)
            :header_(std::move(header)),data_(data){
            /*
            std::transform(data.begin(),data.end(),std::back_inserter(data_),[](std::vector<std::string>& str){
                std::vector<std::string> row;
                for(auto& s:str)
                    row.push_back(boost::trim_copy(s));
                return row;
            });*/

        }




        template<class T>
        bool to(std::vector<std::vector<T>>& num_matrix){

            num_matrix.clear();
            for(auto &vec:data_){
                std::vector<T> num_vec;
                for(auto& s:vec){
                    auto value = T{};
                    try {
                        value = boost::lexical_cast<T>(s);
                    }
                    catch (boost::bad_lexical_cast &e) {
                        std::cout << "Exception caught : " << e.what() << std::endl;
                        return false;
                    }
                    num_vec.push_back(value);
                }
                num_matrix.push_back(num_vec);

            }
            return true;
        };



        bool toDM(casadi::DM& dm){
            std::vector<std::vector<double>> m;
            if(to<double>(m)){
                dm = casadi::DM(m);
                return true;
            }
            return false;
        }


    private:
        std::vector<std::string> header_;
        std::vector<std::vector<std::string>> data_;

    };


    class CSVReader {
    public:
        CSVReader() = delete;
        CSVReader(std::string filename,std::string delimeter=",",bool including_header = false)
            :filename(filename),
             delimeter(delimeter),
            including_header(including_header){}


        CSVData read(){
            std::ifstream file(filename);
            std::vector<std::string> header;
            std::vector<std::vector<std::string> > dataList;
            std::string line = "";

            if(including_header){
                getline(file, line);
                line.erase(line.find_last_not_of(" \n\r\t")+1);
                boost::algorithm::split(header, line, boost::is_any_of(delimeter));
            }

            // Iterate through each line and split the content using delimeter
            while (getline(file, line))
            {
                std::vector<std::string> vec;
                line.erase(line.find_last_not_of(" \n\r\t")+1);
                boost::algorithm::split(vec, line, boost::is_any_of(delimeter));
                for(auto& s:vec){
                    boost::trim(s);
                }
                dataList.push_back(vec);
            }
            // Close the File
            file.close();
            return CSVData{header,dataList};
        }
    private:
        std::string filename;
        std::string delimeter;
        bool including_header;
    };

}

#endif //RACECARPLANNER_CSV_READER_HPP
