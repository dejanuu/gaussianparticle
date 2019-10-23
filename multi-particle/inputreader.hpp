#ifndef INPUTREADER_H
#define INPUTREADER_H

#include "definitions.hpp"
#include <string>
#include <memory>
#include <fstream>
#include <sstream>


class InputReader{
    public:
        InputReader(std::string fname){
            std::string line;
            std::ifstream input(fname.c_str());
            while(std::getline(input,line)){
                int pos = line.find(" "); 
                std::string a = line.substr(0,pos);     
                std::string b = line.substr(pos); 
                mapping.insert(std::pair<std::string,std::string>(a,b));
            }
            input.close();
        }

        ~InputReader(){}

        template <typename T> 
        const T extract (std::string str,const T defaultValue) const {
            T num;
            auto it = mapping.find(str);

            

            if(it == mapping.end()){
                std::cout << str << " was not given, using default value: " << defaultValue << std::endl;
                return defaultValue;
            }else{
                std::istringstream ss(it->second);
                ss >> num;
                std::cout << str << " = " << num << std::endl;
                return num;
            }
        }

        
        const bool arg_found (std::string str) const {
            return (mapping.find(str) != mapping.end());
        }


    private:
        std::map<std::string, std::string> mapping;
};



#endif