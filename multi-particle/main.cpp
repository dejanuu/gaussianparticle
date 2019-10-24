#include <iostream>
#include <list>
#include "inputreader.hpp"
#include "meshreader.hpp"
#include "rng.hpp"
#include "definitions.hpp"
#include "materials.hpp"
#include "physicsengine.hpp"
#include "detector.hpp"
#include "geometry.hpp"
#include "outputwriter.hpp"
#include "tracer.hpp"

int main(int argc, char* argv[]){
    //std::cout << "General Siris: Git Hash: " << g_GIT_SHA1 << std::endl;
    //Set precision of the output
    std::cout << std::setprecision(17);
    if(argc<=1){
        std::cout << "Input file was not given" << std::endl;
        return EXIT_FAILURE;
    }


    std::cout << create_heading(std::string("READING INPUTS FROM: ")+std::string(argv[1]),HEADINGS::HEADING) << std::endl;

    //Create input reader
    InputReader reader(argc,argv);
    
    //Read seed from the input file, -1 means that seed will be generated
    RNG rng(reader.extract<int>("seed",-1));

    std::vector<std::tuple<int,int>> material_inds;
    std::vector<Point> points;
    std::list<Triangle> triangles;


    if(!MeshReader::read_file(reader.extract<std::string>("mesh",std::string("test.off")),
                            points,triangles,material_inds,
                            reader.extract<bool>("centralize",false))) return EXIT_FAILURE;
    

    Materials materials(reader);

    if(!materials.connect_material(material_inds, triangles)) return EXIT_FAILURE;

    Detector detector(reader,rng);
    Geometry geometry(triangles);
    MeshStats meshStats(points);
    PhysicsEngine physEngine(reader, rng, meshStats);



    std::cout << create_heading("TRACER: INITIALIZE",HEADINGS::HEADING) << std::endl;
    Tracer tracer(reader);
    std::cout << create_heading("TRACER: START",HEADINGS::HEADING) << std::endl;
    try{
        tracer.start(geometry, physEngine, materials, detector); 
    }catch(std::exception& e){
        std::cout << create_heading("#!##!#TRACER: ERROR#!#!#!#!#",HEADINGS::HEADING) << std::endl;
        std::cout << e.what() << std::endl;
    }
    std::cout << create_heading("TRACING: FINISHED",HEADINGS::HEADING) << std::endl;
 

    OutputWriter::write_output(meshStats,reader,detector,rng.get_seed());
    std::cout << create_heading("EXIT",HEADINGS::HEADING) << std::endl;
    return EXIT_SUCCESS;
}
