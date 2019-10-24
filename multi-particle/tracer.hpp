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
#include "boost/date_time/posix_time/posix_time_types.hpp"


class Tracer{
    public:
        Tracer(const InputReader& inputreader);

        ~Tracer(){};
        Tracer(const Tracer& other){};
        Tracer& operator=(const Tracer& other){return *this;};
        
        //void process_pass(FortranInterface* fortranInterface,Detector* detector, int nrays);
        void start(Geometry& geometry, PhysicsEngine& physEngine, Materials& materials, Detector& detector); 
        bool allow_continue(int max_rays, int traced_rays,const boost::posix_time::ptime& start_time, int check_interval, double allocated_time, Detector& detector);
    private:
        int nrays = 10;
        int check_time = 10000;
        double allocated_time = 99999.0;

};

#endif