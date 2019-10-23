#ifndef OFFREADER_H
#define OFFREADER_H

#include "definitions.hpp"
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <fstream>
#include <string>
#include <vector>
#include <tuple>

class OFFReader{
    public:
        static void read_file(std::string& fname, std::vector<Point>& points, std::list<Triangle>& triangles, std::vector<std::tuple<int, int>>& material_inds){
            std::ifstream input(fname.c_str());
            if (!input) throw std::invalid_argument(std::string("Cannot open file: ")+std::string(fname));

            std::string tmp;
            int nvert, nfaces, something;
            
            std::getline (input,tmp);
            std::getline (input,tmp);
            while(tmp.find('#') != std::string::npos){
                std::getline (input,tmp);
            }
            while(tmp.size()==0){
                std::getline (input,tmp);
            }
            std::istringstream buffer(tmp);
           
            buffer >> nvert >> nfaces >> something;
            std::cout << "n_vertices: " << nvert << ", n_faces: " << nfaces << std::endl;
            double c1,c2,c3;
            for(int i=0;i<nvert;++i){
                input >> c1 >> c2 >> c3;
                points.push_back(Point(c1,c2,c3));
            }
            
            centralize(points);
            input.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            std::vector<int> face;
            for(int i=0;i<nfaces;++i){
                std::getline(input,tmp);
                std::stringstream ss(tmp);
                std::string item;
                char delim = ' ';
                while (!ss.eof()) { 
                    ss >> item; 
                    face.push_back(std::stoi(item));
                } 

                //while(std::getline(ss, item, delim)) {
                //    std::cout << item << std::endl;
                //    face.push_back(std::stoi(item));
                //}
                
                if(face[0]!=3) throw std::invalid_argument("ERROR while reading an off file: The mesh needs to be triangulized");
                if(face.size()>4){
                    material_inds.push_back(std::make_tuple(face[4],face[5]));
                }else{
                    material_inds.push_back(std::make_tuple(0,1));
                }
                triangles.push_back(Triangle(points[face[1]],points[face[2]],points[face[3]]));
                face.clear();
            }
            input.close();
        }

};


#endif
