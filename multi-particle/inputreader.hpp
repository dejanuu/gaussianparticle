#ifndef INPUTREADER_H
#define INPUTREADER_H

#include "definitions.hpp"
#include <string>
#include <memory>
#include <fstream>
#include <sstream>
#include <algorithm>

class InputReader{
    public:
        InputReader(int argc, char* argv[]){
            std::string line;
            std::ifstream input(argv[1]);
            while(std::getline(input,line)){
                int pos = line.find(" "); 
                std::string a = line.substr(0,pos);     
                std::string b = line.substr(pos); 
                mapping.insert(std::pair<std::string,std::string>(a,b));
            }

            if(argc>2){
                std::cout << "overriding input parameters from terminal" << std::endl;
                for(int i=2;i<argc;i=i+2){
                    std::string parameter = std::string(argv[i]);
                    parameter.erase(std::remove(parameter.begin(), parameter.end(), '-'), parameter.end());
                    std::string value = std::string(argv[i+1]);
                    std::cout << parameter << " " << value << std::endl;
                    auto it = mapping.find(parameter); 
                    if(it == mapping.end()){
                        mapping.insert(std::pair<std::string,std::string>(parameter,value));
                    } else {
                        it->second = value;
                    }
                }
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