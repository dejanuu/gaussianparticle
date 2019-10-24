#ifndef PHYSICSENGINE_H
#define PHYSICSENGINE_H


#include "sray.hpp"
#include "definitions.hpp"
#include "rng.hpp"
#include "misc.hpp"
#include "detector.hpp"
#include "inputreader.hpp"
#include "materials.hpp"
#include "smaterial.hpp"
#include <cmath>
#include <exception>


extern"C" {
void scatter_surface_(double* KE1, double* KF1, double* HL1, double* HR1,
                double* MA1, double* F1, int* preserve1, double* N, double* m1re, double* m1im, 
                double* m2re,double* m2im, double* KE2, double* KF2, double* HL2, double* HR2,
                double* MA2, double* F2, int* preserve2, int* prevent_TR);
}



struct PhysicsEngine{
    private:
        double _max_dist;
        RNG& rng;
        double I_cutoff_limit = 0.001;
        int killswitch_start = 100;
        int max_scattering = 10;
        double wavelen = 2*M_PI;
        double mesh_scale = 1.0;
        double beam_radius = -1.0;
        bool ignore_absorption = false;
        int prevent_TR = 0;


    public:
        PhysicsEngine(InputReader& reader, RNG& _rng, MeshStats& meshStats) : rng(_rng){
            _max_dist = meshStats.max_dist;
            beam_radius = _max_dist;

            ignore_absorption = reader.extract<bool>("ignore_absorption",ignore_absorption);
            prevent_TR = reader.extract<int>("prevent_TR",prevent_TR);

            max_scattering = reader.extract<double>("max_scattering",max_scattering);
            I_cutoff_limit = reader.extract<double>("I_cutoff_limit",I_cutoff_limit);
            killswitch_start = reader.extract<double>("killswitch_start",killswitch_start);
            wavelen = reader.extract<double>("wavelen",wavelen);
            mesh_scale = reader.extract<double>("mesh_scale",mesh_scale);
            double tmp = reader.extract<double>("beam_radius",-1.0);
            if(tmp>=0.0){
                std::cout << "beam_radius > 0: custom beam size is in use" << std::endl;
                beam_radius = tmp/mesh_scale;
            }
            meshStats.compute_area(beam_radius*mesh_scale); 
        }


    

        void scatter_surface(
                            SRay& ray1, SRay& ray2, double* N, Material* material1, Material* material2, 
                            int& preserve1, int& preserve2, int& prevent_TR){

            
            ray2.material = (ray1.material==material1 ? material2 : material1);
            scatter_surface_(ray1.K, ray1.KF, ray1.HL, ray1.HR,
                ray1.M, ray1.F, &preserve1, N, &ray1.material->real, &ray1.material->imag, 
                &ray2.material->real, &ray2.material->imag,  ray2.K, ray2.KF, ray2.HL, ray2.HR,
                ray2.M, ray2.F, &preserve2, &prevent_TR);
        }


        void diffuse_scattering(std::vector<SRay>& stack_diffuse, Detector& detector){
            SRay& ray = stack_diffuse.back();
            ray.property=RAYPROPERTY::DIFFUSESCATTERING;
            detector.register_diffuse_absorption(ray.material->diffuse.diffuse_absorption(ray.F));
            ++ray.n_scattering;
            std::vector<double> rands(2);
            rng.rand(rands,2);
            ray.material->diffuse.scatter(ray.F, ray.K, ray.KF, ray.HL, ray.HR, &rands[0]);

            if(kill_ray(1,ray,max_scattering,I_cutoff_limit,detector)){
                stack_diffuse.pop_back();
            }
        }


        bool escape_diffuse(
                            std::vector<SRay>& stack,std::vector<SRay>& stack_diffuse, 
                            Material* material1, Material* material2, const Point* p, 
                            double* N, Detector& detector){

            SRay& ray1 = stack_diffuse.back();
            double len = generate_dist(ray1);
            double len2=0.0;
            for(int i = 0; i<3;++i) len2 += ((ray1.X[i]-p->cartesian(i))*(ray1.X[i]-p->cartesian(i)));
            if(len2<len*len){
                handle_nondiffuse_scattering(stack_diffuse, stack, stack_diffuse, material1, material2, p , N, detector);
                return true;
            }else{
                detector.register_absorption(traverse(ray1, len*len, mesh_scale, ignore_absorption));
                if(kill_ray(1,ray1,max_scattering, I_cutoff_limit,detector)){
                    stack_diffuse.pop_back();
                    return true;
                }
                ray1.X[0] += ray1.K[0]*len; 
                ray1.X[1] += ray1.K[1]*len;
                ray1.X[2] += ray1.K[2]*len;
  
                return false;
            }
        }


        
        double traverse(SRay& ray1, double len2, double mesh_scale, bool ignore_absorption){
            double abscf = ray1.material->abscf/wavelen;
            //std::cout << "wee " << abscf << std::endl;
            double qabs = 0.0;
            if(abscf<0.0){
                double kekf=0.0;
                for(int i = 0; i<3;++i){
                    kekf += ray1.K[i]*ray1.KF[i];
                }
                
                double attcf = 0.0;
                double lenabs = std::sqrt(len2)*kekf*mesh_scale;
                
                if(abs(abscf*lenabs)<20.0){
                    attcf = exp(abscf*lenabs);
                }

                if(ignore_absorption){
                    attcf = 1.0;
                }

                //std::cout << "kekf" << kekf << ","<< len2 <<  "," << lenabs << "," << (1.0-attcf) <<  "," << std::endl;
                //std::cout << qabs <<"," << (1.0-attcf) <<"," << ray1.F[0] << std::endl;

                qabs = (1.0-attcf)*ray1.F[0];
                multiply_mat_coeff(ray1.F,attcf,16);
            }
            return qabs;
        }


        inline double generate_dist(SRay& ray){
            return ray.material->generate_propagation_distance(rng.rand())/mesh_scale;
        }



        void handle_nondiffuse_scattering(std::vector<SRay>& current_stack, std::vector<SRay>& stack, std::vector<SRay>& stack_diffuse, 
                                        Material* material1, Material* material2, const Point* p , double* N,
                                        Detector& detector){
            SRay& ray1 = current_stack.back();
            
            ray1.property = RAYPROPERTY::SCATTERED;
            //std::cout << "in" << ray1.F[0] << std::endl;
            double len2 = 0.0;
            for(int i = 0; i<3;++i) len2 += ((ray1.X[i]-p->cartesian(i))*(ray1.X[i]-p->cartesian(i)));
            for(int i = 0; i<3;++i) ray1.X[i] = p->cartesian(i);
            detector.register_absorption(traverse(ray1,len2,mesh_scale,ignore_absorption));
           // std::cout << ray1.F[0] << std::endl;
            if(ray1.F[0]!=ray1.F[0]) throw std::logic_error("error found");
            if(kill_ray(1,ray1,max_scattering, I_cutoff_limit,detector)){
                current_stack.pop_back();
                return;
            }

            ++ray1.n_scattering;
            SRay ray2(ray1,true);

            int preserve1 = 0;
            int preserve2 = 0;
            scatter_surface(ray1, ray2, N, material1,material2,preserve1,preserve2,prevent_TR);


            if(preserve2==3){
                detector.register_surface(ray2.F[0]);
                preserve2=0;
            }

            //if(preserve1) std::cout << "1:" << ray1.F[0] << std::endl;
            //if(preserve2) std::cout << "2:"  <<ray2.F[0] << std::endl;
            if(ray1.n_scattering<killswitch_start){
                if(kill_ray(preserve1,ray1,max_scattering, I_cutoff_limit,detector)){
                    current_stack.pop_back();
                }
               
                if(!kill_ray(preserve2,ray2,max_scattering, I_cutoff_limit,detector)){
                    if(ray2.material->diffuse_inclusions){
                        stack_diffuse.push_back(ray2);
                    }else{
                        stack.push_back(ray2);
                    }
                }
            }else if(preserve1+preserve2==2){
                //std::cout << "ee" << std::endl;
                double F = ray1.F[0]+ray2.F[0];
                if(F*rng.rand()>ray1.F[0]){
                    current_stack.pop_back();
                    multiply_mat_coeff(ray2.F,F/ray2.F[0],16);
                    if(ray2.material->diffuse_inclusions){
                        stack_diffuse.push_back(ray2);
                    }else{
                        stack.push_back(ray2);
                    }
                }else{
                    multiply_mat_coeff(ray1.F,F/ray1.F[0],16);
                }
            }else if(preserve1==0){
                current_stack.pop_back();
                if(ray2.material->diffuse_inclusions){
                    stack_diffuse.push_back(ray2);
                }else{
                    stack.push_back(ray2);
                }   
                
            }
        }

        bool kill_ray(int flag, SRay& ray, const int max_scattering, const double I_cutoff, Detector& detector){
            if(flag==0) return true;
            if(ray.F[0]<I_cutoff || ray.n_scattering>max_scattering){
                detector.register_cutoff(ray.F[0]); 
                return true;
            }
            return false;
        }


        void handle_escaping(std::vector<SRay>& stack, Detector& detector, Material* material_surrouding){
            SRay& ray = stack.back();
            if(ray.property ==RAYPROPERTY::UNSCATTERED){
                detector.register_unscattered();
            }else{
                if(ray.material!=material_surrouding) throw std::logic_error("incorrect material");
                detector.register_scattered(ray.F,ray.K,ray.HL,ray.theta,ray.phi);
            }
            stack.pop_back();
        }

        void generate_rays(std::vector<SRay>& stack, Detector& detector,Materials& materials){
            std::vector<double> rands(4);
            rng.rand(rands,4);
            for(int i=0;i<1;++i){
                double theta = acos(1.0-rands[0]*2.0);
                double phi = M_PI*2.0*rands[1];
            
                double X[3];
                double tmp = rands[2]*2*M_PI;
                double r = sqrt(rands[3])*beam_radius;
                X[0] = r*cos(tmp);
                X[1] = r*sin(tmp);
                X[2] = -_max_dist;

                local_to_world_space(theta,phi,X,0,1,2);
            
                stack.push_back(SRay(wavelen,theta,phi,X,materials.get_surrounding_material()));
  
                detector.register_ray();
            }
        }
};





#endif
