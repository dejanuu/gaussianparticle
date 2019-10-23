#ifndef RAY_H
#define RAY_H


#include <cmath>
#include <stdexcept>
#include "definitions.hpp"
#include "smaterial.hpp"

enum RAYPROPERTY{ UNSCATTERED, SCATTERED, HIT, MISS,DIFFUSESCATTERING};


template <typename T>
struct MRay{
    double X[3];
    double K[3];
    double theta;
    double phi;
    RAYPROPERTY property = UNSCATTERED;
    int n_scattering = 0;
    T primitive;
    Material* material;

    MRay(double theta, double phi, const double* X, Material* material){
        K[0]=sin(theta)*cos(phi);
        K[1]=sin(theta)*sin(phi);
        K[2]=cos(theta);

        this->X[0] = X[0];
        this->X[1] = X[1];
        this->X[2] = X[2];
        this->property = property;
        this->primitive = primitive;

        this->theta = theta;
        this->phi = phi;
        this->material = material;
    } 

    MRay(const MRay& ray){
        std::copy(ray.X,  ray.X + 3, X);
        std::copy(ray.K, ray.K + 3, K);
        n_scattering = ray.n_scattering;
        property = ray.property;
        primitive = ray.primitive;
        theta = ray.theta;
        phi = ray.phi;
        material = ray.material;
    }


    /*
    inline Material* get_material(Material* material1, Material* material2){
        if(material1==material){
            return material2;
        }else if(material2==material){
            return material1;
        }else{
            throw std::domain_error(std::string("Nonmanifold geometry. Check materials and the mesh"));
        }
    }
    */
    inline Ray extract_ray(){
        return Ray(Point(X[0],X[1],X[2]),Vector3(K[0],K[1],K[2]));
    }


};





#endif