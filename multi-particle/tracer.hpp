#ifndef TRACER_H
#define TRACER_H


#include <string>
#include <fstream>
#include <vector>
#include <iterator>
#include "definitions.hpp"
#include "rng.hpp"
#include "inputreader.hpp"
#include "detector.hpp"
#include "geometry.hpp"
#include <math.h>



class Tracer{
    public:
        Tracer(const InputReader& inputreader);

        ~Tracer(){};
        Tracer(const Tracer& other){};
        Tracer& operator=(const Tracer& other){return *this;};
        
        //void process_pass(FortranInterface* fortranInterface,Detector* detector, int nrays);
        void start(Geometry& geometry, PhysicsEngine& physEngine, Materials& materials, Detector& detector); 

    private:
        int nrays = 10;

};

#endif