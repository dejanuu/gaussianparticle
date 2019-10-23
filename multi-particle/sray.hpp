#ifndef SRAY_H
#define SRAY_H

#include "definitions.hpp"
#include "mray.hpp"
#include "misc.hpp"
#include <cmath>
#include <stdexcept>

struct SRay : public MRay<PrimitiveType>{
    double KF[3];
    double F[16];
    double HL[6];
    double HR[6];
    double M[2];
    double wavelen;

    //SRay(double theta, double phi, double* X, Material* material) : MRay<PrimitiveType>(theta,phi,X,material){};

    SRay(double wavelen, double theta, double phi, const double* X, Material* material) :  MRay<PrimitiveType>(theta,phi,X,material) { 
        fill_zero(F,16);
        F[0] = 1.0;
        F[5] = 1.0;
        F[10] = 1.0;
        F[15] = 1.0;
        fill_zero(HL,6);
        fill_zero(HR,6);
        HL[0] = 1.0;
        HR[2] = 1.0;
        std::copy(K, K + 3, KF);
        local_to_world_space(theta,phi,HL,0,2,4);
        local_to_world_space(theta,phi,HR,0,2,4);
        this->wavelen = wavelen;
  
        double mre=material->real;
        double mim=0.0;//material->imag;
        
        //double d = mre*mre-mim*mim;
        //M[0] = std::sqrt(0.5*(d+sqrt(d*d+4.0*std::pow((mre*mim),2.0))));
        //M[1] = std::sqrt(0.5*(-d+sqrt(d*d+4.0*std::pow((mre*mim),2.0))));
        M[0] = mre;
        M[1] = 0.0;
    } 

    SRay(const SRay& ray,bool clean) : MRay<PrimitiveType>(ray){
        if(!clean){
            std::copy(ray.KF, ray.KF + 3, KF);
            std::copy(ray.F,  ray.F + 16, F);
            std::copy(ray.HL, ray.HL + 6, HL);
            std::copy(ray.HR, ray.HR + 6, HR);
            std::copy(ray.M, ray.M + 2, M);
        }  
        this->wavelen = ray.wavelen;
    }
};

#endif