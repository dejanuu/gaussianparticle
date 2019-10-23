#include "tracer.hpp"

Tracer::Tracer(const InputReader& inputReader){
    nrays = inputReader.extract<int>("nrays",nrays);
}

void Tracer::start(Geometry& geometry, PhysicsEngine& physEngine, Materials& materials, Detector& detector){

    std::vector<SRay> stack_diffuse;
    std::vector<SRay> stack;
    stack_diffuse.reserve(100);
    stack.reserve(100);
    int stack_limit = 50;
    while(detector.nrays<nrays){
        physEngine.generate_rays(stack,detector,materials);
        //std::cout << "start new" << std::endl;
        while(stack.size()>0 || stack_diffuse.size()>0){
            while(stack.size()>0){
                //std::cout << "nondiffuse " << stack.size() << " " << stack.back().F[0] << "," << stack.back().X[0]<< "," << stack.back().X[1]<< "," << stack.back().X[2]<< std::endl;
                geometry.trace_in_nondiffuse(stack,stack_diffuse, physEngine, detector, materials);
            }
            while(stack_diffuse.size()>0){
            //    std::cout << "diffuse " << stack_diffuse.size() << std::endl;
                geometry.trace_in_diffuse(stack,stack_diffuse,physEngine, detector, materials);
            }
        }
    }

}

