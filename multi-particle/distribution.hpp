#ifndef DISTRIBUTION_H
#define DISTRIBUTION_H

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <iterator>
#include <sstream> 
#include <cmath>

class Distribution{
    public:
        std::vector<double> distX;
        std::vector<double> cdf;
        std::string fname = std::string();
        
        Distribution() : distX(), cdf(){};
        

        Distribution(const Distribution& o){
            //distX = o.distX;
            //cdf = o.cdf;
            //fname = o.fname;
            distX = std::move(o.distX);
            cdf = std::move(o.cdf);
            fname = o.fname;
        }

        Distribution(Distribution&& o) noexcept{    
            distX = std::move(o.distX);
            cdf = std::move(o.cdf);
            fname = o.fname;
        }
        

        Distribution& operator=(const Distribution& o)
        {
            //std::cout << o.distX.size() << " " << o.cdf.size() << std::endl;
            *this = Distribution(o);
            return *this;
        }    
        

       Distribution& operator=(Distribution&& o) noexcept
        {
            distX = std::move(o.distX);
            cdf = std::move(o.cdf);
            fname = o.fname;
            return *this;
        }

        

        Distribution(std::string _fname) : distX(), cdf(){
            std::ifstream file;
            fname = std::string(_fname);
            std::string line;
            distX.push_back(0);
            cdf.push_back(0);
            file.open(fname.c_str());
            if(!file.is_open()) std::cout << "File does not exist" << std::endl;
            while(getline(file,line) ){
                std::stringstream ss(line);
                double arg1,arg2;
                ss >> arg1 >> arg2;
                distX.push_back(arg1);
                cdf.push_back(arg2); 
            }
            distX.push_back(std::pow(10,7));
            cdf.push_back(1);
            file.close();   

        }

        double pick_dist(double u){
            auto iter = std::lower_bound(cdf.begin(),cdf.end(),u);
            int indx = iter-cdf.begin()-1;
            double x = (u-cdf[indx])/(cdf[indx+1]-cdf[indx])*(distX[indx+1]-distX[indx])+distX[indx];
            //std::cout << x << " " << u << " " << indx << std::endl;
            //std::cout << cdf[indx] <<" "<< cdf[indx+1] <<" "<< distX[indx] << " "<< distX[indx+1] << std::endl;
            return x;
        }


};




#endif