#include "tracer.hpp"

Tracer::Tracer(const InputReader& inputReader){
    nrays = inputReader.extract<int>("nrays",nrays);
    check_time = inputReader.extract<int>("check_time_after_nrays",check_time); 
    if(check_time>0){
        allocated_time = inputReader.extract<double>("allocated_time_in_hours",allocated_time); 
    }
    
    
}


bool Tracer::allow_continue(int max_rays, int traced_rays,const boost::posix_time::ptime& start_time, int check_interval, double allocated_time, Detector& detector){
    //I'm sorry for this
    if(traced_rays>0 && check_interval>0){
        if(traced_rays%check_interval==0){
            boost::posix_time::ptime current_time = boost::posix_time::second_clock::universal_time();
            if((current_time-start_time).total_seconds()/3600.0>allocated_time){
                detector.raise_time_kill();
                return false;
            }
        }
    }
    return (traced_rays<max_rays);
}




void Tracer::start(Geometry& geometry, PhysicsEngine& physEngine, Materials& materials, Detector& detector){
    detector.set_max_nrays(nrays);
    std::vector<SRay> stack_diffuse;
    std::vector<SRay> stack;
    stack_diffuse.reserve(100);
    stack.reserve(100);
    boost::posix_time::ptime start_time = boost::posix_time::second_clock::universal_time();
    while(allow_continue(nrays,detector.nrays,start_time,check_time,allocated_time,detector)){
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

