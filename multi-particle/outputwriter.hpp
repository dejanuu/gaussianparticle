#ifndef OUTPUTWRITER_H
#define OUTPUTWRITER_H

#include "definitions.hpp"
#include "inputreader.hpp"
#include "detector.hpp"
#include <string>
#include <cmath>


class OutputWriter{
    public:
        static void write_output(MeshStats& meshStats, InputReader& reader,Detector& detector){
            try{

                std::string output_file = reader.extract<std::string>("output_file",std::string("outputS.out"));
                std::string output_file2 = reader.extract<std::string>("pmatrix_out",std::string("pmatrix.out"));
                std::string details = reader.extract<std::string>("details_out",std::string("details.out"));

                double mesh_scale = reader.extract<double>("mesh_scale",1.0);
                double _max_dist = meshStats.max_dist;
                double area = M_PI*std::pow(mesh_scale*_max_dist,2);


                int nbins = detector.nbins;
                int nrays = detector.nrays;


                write_details(nbins,nrays, detector,area);

                std::ofstream out(details.c_str());
                std::streambuf *coutbuf = std::cout.rdbuf();
                std::cout.rdbuf(out.rdbuf()); 
                write_details(nbins,nrays, detector,area);
                std::cout.rdbuf(coutbuf); 


                std::vector<double>& thetas = detector.thetas;
                double* F = detector.normalize_output();
                             

                std::ofstream myfile;
                myfile.open(output_file.c_str());
                for(int i=0;i<nbins+1;++i){
                    myfile << thetas[i]*180/M_PI << " ";
                    for(int k=0;k<4;++k){
                        for(int j=0;j<4;++j){
                            myfile << F[i*16+k+4*j] << " ";   
                        }
                    }      
                    myfile << std::endl;
                }
                myfile.close();


                myfile.open(output_file2.c_str());
                for(int i=0;i<nbins+1;++i){
                    myfile << thetas[i]*180/M_PI << " ";
                    for(int j=0;j<16;++j){
                        if(j==0 || j==1 || j==5 || j==10 || j==11 || j==15 ){
                            myfile << F[i*16+j] << " ";   
                        }
                    }         
                    myfile << std::endl;
                }
                myfile.close();


            }catch (std::exception e) {
                std::cout << "Output writer: " << e.what() << std::endl;
            }
        }

    private:
        static void write_details(int nbins,int nrays, Detector& detector, double area){
            std::cout << create_heading("OUTPUT",HEADINGS::HEADING) << std::endl;
            std::cout << "From: " << nrays << " rays " << detector.n_unhit << " unhit" << std::endl;
            std::cout << "Area: " << area*(nrays-detector.n_unhit)/nrays << std::endl;
            std::cout << create_heading("",HEADINGS::SUB_HEADING) << std::endl;
            int tmp = nrays-detector.n_unhit;
            std::cout << "Q_abs      : " << detector._qabs/tmp << std::endl;
            std::cout << "Q_abs_diff : " << detector._qabs_diff/tmp << std::endl;
            std::cout << "Q_sca      : " << detector._qsca/tmp << std::endl;
            std::cout << "Q_stop     : " << detector._qstop/tmp << std::endl;
            std::cout << "Q_surf     : " << detector._qsurf/tmp << std::endl;
            std::cout << create_heading("",HEADINGS::SUB_HEADING) << std::endl;
            std::cout << "Total: " << (detector._qstop +detector._qabs+detector._qabs_diff+detector._qsca +detector._qsurf)/tmp << std::endl;
            std::cout << create_heading("SAVING OUTPUT",HEADINGS::SUB_HEADING) << std::endl;
        }

};


#endif