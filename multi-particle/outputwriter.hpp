#ifndef OUTPUTWRITER_H
#define OUTPUTWRITER_H

#include "definitions.hpp"
#include "inputreader.hpp"
#include "detector.hpp"
#include <string>
#include <cmath>


class OutputWriter{
    public:
        static void write_output(MeshStats& meshStats, InputReader& reader,Detector& detector, int seed){
            try{

                std::string output_file = reader.extract<std::string>("output_file",std::string("outputS.out"));
                std::string output_file2 = reader.extract<std::string>("pmatrix_out",std::string("pmatrix.out"));
                std::string details = reader.extract<std::string>("details_out",std::string("details.out"));

               

                double area = meshStats.area;


                int nbins = detector.nbins+detector.nbins_fine_details;
                int nrays = detector.nrays;


                write_details(nbins,nrays, detector,area,seed);

          
                std::ofstream out(details.c_str());
                std::streambuf *coutbuf = std::cout.rdbuf();
                std::cout.rdbuf(out.rdbuf()); 
                write_details(nbins,nrays, detector,area,seed);
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
                        if(j==0 || j==4 || j==5 || j==10 || j==14 || j==15 ){
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
        static void write_details(int nbins,int nrays, Detector& detector, double area, int seed){
            std::cout << create_heading("OUTPUT",HEADINGS::HEADING) << std::endl;
            std::cout << "From: " << nrays << " rays " << detector.n_unhit << " unhit " << detector.direct_transmission << " hit but not scattered" << std::endl;
            std::cout << "Real_Area: " << area*(nrays-detector.n_unhit)/nrays << ", Eff_area: " << area*(nrays-detector.n_unhit-detector.direct_transmission)/nrays << std::endl;
            std::cout << create_heading("",HEADINGS::SUB_HEADING) << std::endl;
            int tmp = nrays-detector.n_unhit-detector.direct_transmission;
            std::cout << "Q_abs      : " << detector._qabs/tmp << std::endl;
            std::cout << "Q_abs_diff : " << detector._qabs_diff/tmp << std::endl;
            std::cout << "Q_sca      : " << detector._qsca/tmp << std::endl;
            std::cout << "Q_stop     : " << detector._qstop/tmp << std::endl;
            std::cout << "Q_surf     : " << detector._qsurf/tmp << std::endl;
            std::cout << create_heading("",HEADINGS::SUB_HEADING) << std::endl;
            std::cout << "Total      : " << (detector._qstop +detector._qabs+detector._qabs_diff+detector._qsca +detector._qsurf)/tmp << std::endl;
            std::cout << create_heading("",HEADINGS::SUB_HEADING) << std::endl;
            std::cout << "Seed       : " << seed << std::endl;
            if(detector.time_kill) std::cout <<            "Time kill  : " << detector.nrays-detector.n_unhit-detector.direct_transmission << " " << detector.nrays_max-detector.n_unhit-detector.direct_transmission << std::endl;
            if(detector.force_interaction==1) std::cout << "Q_dt       : " << detector._qdt/detector.direct_transmission  << std::endl;
            std::cout << create_heading("SAVING OUTPUT",HEADINGS::SUB_HEADING) << std::endl;
            //std::cout << "General Siris: Git Hash: " << g_GIT_SHA1 << std::endl;        
        }

};


#endif
