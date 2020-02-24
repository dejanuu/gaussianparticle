#include "offreader.hpp"

//To do: separate the material reading to make this code more general. Some other function
//can decode the color input for SIRIS
//Try to use std::vector instead list, but requires that the CGAL Tree supports this

//Read OFF file format
void OFFReader::read_file(std::string& fname, std::vector<Point>& vertices, 
                                std::list<Triangle>& triangles, std::vector<std::tuple<int, int>>& material_inds, 
                                bool centralize_vertices){
    //establish input stream
    std::ifstream input(fname.c_str());
    //check whether the file is readable
    if (!input) throw std::invalid_argument(std::string("Cannot open file: ")+std::string(fname));

    std::string tmp;
    int nVerts, nFaces, something;

    //reserve memory (do not spend time allocating memory while reading. Important for larger meshes)
    vertices.reserve(nVerts);

    //Get rid of the "OFF"
    std::getline (input,tmp);
    //Read all the unnecessary information from the start of the file (#) until reaching the line that does not contain #
    std::getline (input,tmp);
    while(tmp.find('#') != std::string::npos){
        std::getline (input,tmp);
    }
    //Get rid of all the empty lines
    while(tmp.size()==0){
        std::getline (input,tmp);
    }
    //Now we should have reached the line that contains nVerts,nFaces,nEdges 
    std::istringstream buffer(tmp);
    //nEdges can be anything
    buffer >> nVerts >> nFaces >> something;
    std::cout << "n_vertices: " << nVerts << ", n_faces: " << nFaces << std::endl;
    double c1,c2,c3;
    //store all vertices
    for(int i=0;i<nVerts;++i){
        input >> c1 >> c2 >> c3;
        vertices.push_back(Point(c1,c2,c3));
    }

    //Centralize the mesh. The second centralization is just for checkin purposes
    if(centralize_vertices){
        centralize(vertices);
        centralize(vertices);
    }

    //ignore the new line symbol as a word 
    input.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    std::vector<int> faces;

    //read vertices
    for(int i=0;i<nFaces;++i){
        std::getline(input,tmp);
        //trim input, makes life so much easier
        boost::trim(tmp);
        std::stringstream ss(tmp);
        std::string item;
        char delim = ' ';
        //split the line to elements
        while (!ss.eof()) { 
            ss >> item;
            faces.push_back(std::stoi(item));
        } 

        //while(std::getline(ss, item, delim)) {
        //    std::cout << item << std::endl;
        //    face.push_back(std::stoi(item));
        //}
        
        //There Siris compatible meshes contain only triangles
        if(faces[0]!=3) throw std::invalid_argument("ERROR while reading an off file: The mesh needs to be triangulized");
        
        //Check whether color coding is present
        if(faces.size()>4){
            //It was
            material_inds.push_back(std::make_tuple(faces[4],faces[5]));
        }else{
            //Set default (1 material)
            material_inds.push_back(std::make_tuple(0,1));
        }
        //Let's store the faces to the mesh's triangles
        triangles.push_back(Triangle(vertices[faces[1]],vertices[faces[2]],vertices[faces[3]]));
        //clear temporary data about the faces
        faces.clear();
    }
    //close the stream
    input.close();
}

