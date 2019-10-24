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
void finalize_ray_(double *F, double* K,  double* HL, double *theta, double* phi, double * theta_low_limit, double * theta_upper_limit, double* rands);
}





class Detector{

    public:
        Detector(const InputReader& inputReader, RNG& _rng) : rng(_rng){
            nbins= inputReader.extract<int>("nbins",nbins);
            nbins_fine_details_start= inputReader.extract<double>("nbins_fine_details_start",nbins_fine_details_start);
           
            force_interaction= inputReader.extract<int>("force_interaction",force_interaction);
            if(nbins_fine_details_start<180.0){
                nbins_fine_details= inputReader.extract<int>("nbins_fine_details",nbins_fine_details);
            }else{
                nbins_fine_details = 0;
            }
            if(nbins_fine_details_start<=90) std::cout << "WARNING: nbins_fine_details_start<90" << std::endl;
            nbins_fine_details_start = nbins_fine_details_start/180*M_PI;
           
            int total_bins = nbins_fine_details+nbins;
            thetas.reserve(total_bins+1);
            F = new double[16*(total_bins+1)];
            fill_zero(F,16*(total_bins+1));

            for(int i=0;i<nbins;++i){
                thetas.push_back(i*(nbins_fine_details_start)/(nbins));
            }
            for(int i=0;i<nbins_fine_details;++i){
                thetas.push_back(nbins_fine_details_start+i*(M_PI-nbins_fine_details_start)/(nbins_fine_details));
            }
            thetas.push_back(M_PI);
            lower_limit = (thetas[0]+thetas[1])*0.5;
            upper_limit = (thetas[total_bins]+thetas[total_bins-1])*0.5;

            if(thetas.size()!=total_bins+1) std::cout << "WARNING: total_bins and n thetas don't match";
        }

        ~Detector(){
            delete[] F;
        };

        void register_scattered(double* F, double* K, double* HL, double theta, double phi){
            

            double rand = rng.rand();
            finalize_ray_(F, K, HL, &theta, &phi, &lower_limit, &upper_limit, &rand);

            //double bin0 = std::cos(0.5 * M_PI/nbins);
            //KOUT(3)<bin0
            //std::cout << "ready" << std::endl;
            int kbin = 0;
            if(K[2]<std::cos(upper_limit)){
                kbin = nbins+nbins_fine_details;
                //std::cout << thetas[kbin] << " " << std::acos(K[2]);
            }else if(K[2]<=std::cos(nbins_fine_details_start)){
                //std::cout << std::acos(K[2]) << std::endl;
                kbin = nbins+ std::floor(0.5+(std::abs(std::acos(K[2])-nbins_fine_details_start))/((M_PI-nbins_fine_details_start)/nbins_fine_details)); 
                //std::cout << thetas[kbin] << " " << std::acos(K[2]) << thetas[kbin+1];
            }else if(K[2]<=std::cos(lower_limit)){
                kbin = std::floor(0.5+std::acos(K[2])/(nbins_fine_details_start/nbins));  
                //std::cout << thetas[kbin] << " " << std::acos(K[2]);
            }

            if(kbin==0 && force_interaction==1 && F[0]>=0.9999999){
                _qdt+=F[0];
                direct_transmission+=1;
                return;
            }

            _qsca += F[0];
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

        inline void raise_time_kill(){
            time_kill=true;
        }

        inline void set_max_nrays(int nrays){
            nrays_max=nrays;
        }

        


        RNG& rng;
        double lower_limit;
        double upper_limit;
        double nbins_fine_details_start = 180;
        int nbins_fine_details = 0;
        double _qabs = 0.0;
        double _qabs_diff = 0.0;
        double _qsca = 0.0;
        double _qstop = 0.0;
        double _qsurf = 0.0;
        double _qdt = 0.0;
        bool time_kill = false;
        int direct_transmission=0;
        int n_unhit=0;
        int nbins=180;
        int force_interaction = 0;
        int nrays=0;
        int nrays_max=0;
        double * F;
        std::vector<double> thetas;


        double* normalize_output(){
            double dTheta =  lower_limit;
            double tmp = 2.0/(_qsca/(1.0*nrays-n_unhit-direct_transmission));
            double norm = tmp/((1.0*nrays-n_unhit-direct_transmission)*(1.0-cos(dTheta)));
            multiply_mat_coeff(F, norm, 16);

            dTheta = upper_limit;
            norm =  tmp/((1.0*nrays-n_unhit-direct_transmission)*(cos(dTheta)+1.0));
            multiply_mat_coeff(&F[16*(nbins+nbins_fine_details)], norm, 16);


            

            for(int i=1;i<nbins+nbins_fine_details;++i){
                double theta1 =  0.5*(thetas[i]+thetas[i-1]);
                double theta2 =  0.5*(thetas[i+1]+thetas[i]);
                norm = tmp/((1.0*nrays-n_unhit-direct_transmission)*(cos(theta1)-cos(theta2)));
                //norm = 1.0/((2.0*M_PI/nbins)*(cos(0.5*(thetas[i]+thetas[i-1]))-cos(0.5*(thetas[i+1]+thetas[i]))));
                multiply_mat_coeff(&F[i*16], norm, 16);
                //multiply_mat_coeff(&F[i*16], 1.0, 16);
            } 

            return F;
        }

};



#endif