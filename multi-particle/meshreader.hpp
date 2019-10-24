#ifndef MESHREADER_H
#define MESHREADER_H

#include "definitions.hpp"
#include "offreader.hpp"
#include <cmath>
#include <sstream>
#include <iostream>
#include <tuple>

//Stack overflow kek, but I want to go home. it's saturday morning 00:27
inline bool hasEnding(std::string const &fullString, std::string const &ending) {
    if (fullString.length() >= ending.length()) {
        return (0 == fullString.compare (fullString.length() - ending.length(), ending.length(), ending));
    } else {
        return false;
    }
}





class MeshReader{
    public:
        static bool read_file(std::string fname, std::vector<Point>& points, std::list<Triangle>& triangles, std::vector<std::tuple<int, int>>& material_inds, bool centralize){
            std::cout << create_heading(std::string("READ MESH FROM: ")+fname,HEADINGS::HEADING) << std::endl;
            try{
                if(hasEnding(fname,".off")){
                    OFFReader::read_file(fname,points,triangles,material_inds,centralize);
                }else{
                    throw std::invalid_argument("File format not supported"); 
                } 
                std::cout << create_heading("MESH READ SUCCESFULLY",HEADINGS::SUB_HEADING) << std::endl;
            }catch (const std::exception& e) { 
                std::cout << e.what() << std::endl; 
                return false;
            }

            return true;
        }

    void init_periodic(std::vector<Point>& points){
        std::vector<Point>::iterator it; 
        double maxv = -MAXDBLVALUE;
        for(it = points.begin(); it != points.end(); it++)    {
            if((*it)[2]>maxv) maxv = (*it)[2];
        }
        
        for(it = points.begin(); it != points.end(); it++)    {
            Point & point = *it;
            point = point-Vector3(0.0,0.0,maxv);;
        }
    }





};

#endif