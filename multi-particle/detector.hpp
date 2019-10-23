#ifndef DETECTOR_H
#define DETECTOR_H

#include "definitions.hpp"
#include "misc.hpp"
#include "inputreader.hpp"
#include "sray.hpp"
#include "mray.hpp"
#include "rng.hpp"
#include <cmath>

extern"C" {
void finalize_ray_(double *F, double* K,  double* HL, double *theta, double* phi, int* nbins, double* rands);
}





class Detector{

    public:
        Detector(const InputReader& inputReader, RNG& _rng) : rng(_rng){
            nbins= inputReader.extract<int>("nbins",nbins);
            thetas.reserve(nbins+1);
            F = new double[16*(nbins+1)];
            fill_zero(F,16*(nbins+1));
            for(int i=0;i<nbins+1;++i){
                thetas.push_back(i*M_PI/(nbins));
            }
        }

        ~Detector(){
            delete[] F;
        };

        void register_scattered(double* F, double* K, double* HL, double theta, double phi){
            _qsca += F[0];

            double rand = rng.rand();
            finalize_ray_(F, K, HL, &theta, &phi, &nbins, &rand);

            double bin0 = std::cos(0.5 * M_PI/nbins);
            //KOUT(3)<bin0
            int kbin = 0;
            if(K[2]<-bin0){
                kbin = nbins;
            }else if(K[2]<=bin0){
                kbin = std::floor(0.5+std::acos(K[2])/(M_PI/nbins));  
            }
            add_to_mat(this->F,F,kbin,16,16);
        }

        inline void register_diffuse_absorption(double qabs_diff){
            _qabs_diff += qabs_diff;
        }

        inline void register_absorption(double qabs){
            _qabs+=qabs;
        }

         inline void register_cutoff(double qstop){
            _qstop+=qstop;
        }

        inline void register_unscattered(){
            ++n_unhit;
        }

        inline void register_ray(){
            ++nrays;
        }

        inline void register_surface(double qsurf){
            _qsurf +=qsurf;
        }

        RNG& rng;
        double _qabs = 0.0;
        double _qabs_diff = 0.0;
        double _qsca = 0.0;
        double _qstop = 0.0;
        double _qsurf = 0.0;
        int n_unhit=0;
        int nbins=180;
        int nrays=0;
        double * F;
        std::vector<double> thetas;

        double* normalize_output(){
            double bin = M_PI/(1.0*nbins);
            double norm = 1.0/((1.0*nrays-n_unhit)*(1.0-cos(0.5*bin)));
            multiply_mat_coeff(F, norm, 16);

            norm =  1.0/((1.0*nrays-n_unhit)*(cos((nbins-0.5)*bin)+1.0));
            multiply_mat_coeff(&F[16*nbins], norm, 16);


            

            for(int i=1;i<nbins;++i){
                
                norm = 1.0/((1.0*nrays-n_unhit)*(cos((i-0.5)*bin)-cos((i+0.5)*bin)));
                //norm = 1.0/((2.0*M_PI/nbins)*(cos(0.5*(thetas[i]+thetas[i-1]))-cos(0.5*(thetas[i+1]+thetas[i]))));
                multiply_mat_coeff(&F[i*16], norm, 16);
                //multiply_mat_coeff(&F[i*16], 1.0, 16);
            } 

            return F;
        }

};



#endif