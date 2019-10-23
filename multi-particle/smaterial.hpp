
#ifndef SMATERIAL_H
#define SMATERIAL_H


#include "inputreader.hpp"
#include "definitions.hpp"
#include "misc.hpp"

#include <sstream>
#include <iostream>
#include <stdexcept>
#include <fstream>
#include <map>
#include <tuple>


extern"C" {
void create_diffuse_scatterer_(double* P, double* CSRN, double* XP, double* coeffs, int* np, int* nrn);
};

extern"C" {
void scatter_diffuse_(double* F, double* K1, double* KF1, double* HL, double* HR,
                    double* XP, double* CSRN, double* coeffs, int* np, int* nrn, double* rands);
};




struct DiffuseMaterial
{
    double _mean_free_path = 1.0;
    double _albedo = 1.0;
    std::string _fname;
    int tmp_indx;
    int np;
    double* data = nullptr;
    double* csrn= nullptr;
    double* xp = nullptr;
    double* coeffs = nullptr;
    int nrn = 100;
    int reserved_block = (nrn+1+np+np*6);


    DiffuseMaterial(){
    }

    DiffuseMaterial(double mfp, double albedo, std::string fname) : _fname(fname){
        _mean_free_path = mfp;
        _albedo = albedo;
    }

    DiffuseMaterial(const DiffuseMaterial &o) {
        _mean_free_path = o._mean_free_path;
        _albedo = o._albedo;
        _fname = o._fname;
        tmp_indx = o.tmp_indx;
        np = o.np;
        data = o.data;
        csrn = o.csrn;
        xp = o.xp;
        coeffs = o.coeffs;
        nrn = o.nrn;
        reserved_block = o.reserved_block; 
    }


    ~DiffuseMaterial(){
    }

    int request_n_entries(std::vector<double>& tmp){
        //Open the File
        std::ifstream in(_fname.c_str());
        if(!in) throw std::invalid_argument(std::string("Couldn't open the file:")+std::string(_fname.c_str()));
        tmp_indx = tmp.size();
        np=0;
        std::string str;
        while(std::getline(in, str)){
            ++np;
            std::stringstream ss(str);
            for(int i=0;i<7;++i){
                double value;
                ss >> value;
                //std::cout << value << " ";
                tmp.push_back(value);
            }
            //std::cout << std::endl;
        }
        in.close();
        reserved_block = (nrn+1+np+4*np*6);
        return reserved_block;
    }

    int initialize_material(std::vector<double>& tmp, double* storage){
        data = storage;
        csrn = storage;
        xp =  &storage[nrn+1];
        coeffs = &storage[nrn+1+np];
        create_diffuse_scatterer_(&tmp[0], csrn, xp, coeffs, &np, &nrn);
        return reserved_block;
    }

    double diffuse_absorption(double* F){
        double qabs = (1.0-_albedo)*F[0];
        multiply_mat_coeff(F, _albedo, 16);
        return qabs;
    }

    void scatter(double * F, double* K, double* KF, double* HL, double* HR, double* rands){
        scatter_diffuse_(F,K,KF,HL,HR,xp,csrn,coeffs,&np,&nrn,rands);
    }

};


struct Material{
    public:
        double real = 1.0;
        double imag = 0.0;
        double abscf = 0.0;
        bool diffuse_inclusions = false;
        DiffuseMaterial diffuse;

        Material(){
        }

        Material(double _real, double _imag){
            check_arguments(_real,_imag);
            real = _real;
            imag = _imag;
        }

        Material(const Material& o){
            real = o.real;
            imag = o.imag;
            diffuse_inclusions = o.diffuse_inclusions;
            abscf = o.abscf;
            if(diffuse_inclusions) diffuse = o.diffuse;
        }

        Material(Material&& o) noexcept{
            real = o.real;
            imag = o.imag;
            diffuse_inclusions = o.diffuse_inclusions;
            abscf = o.abscf;
            if(diffuse_inclusions) diffuse = o.diffuse;
        }


        Material& operator=(const Material& o)
        {
            return *this = Material(o);
        }

        Material& operator=(Material&& o) noexcept
        {
            real =  o.real;
            imag =  o.imag;
            diffuse_inclusions = o.diffuse_inclusions;
            if(diffuse_inclusions) diffuse = o.diffuse;
            return *this;
        }

        ~Material(){
        }

        friend std::istream& operator >> (std::istream & input, Material& o){
            input >> o.real  >> o.imag >> o.diffuse_inclusions;
            o.check_arguments( o.real,o.imag);
            if(o.diffuse_inclusions){
                double mean_free_path,albedo;
                std::string fname;
                input >> mean_free_path >> albedo >> fname;
                o.diffuse = DiffuseMaterial(mean_free_path,albedo,fname);
            }
            
        }

        friend std::ostream& operator<<(std::ostream & stream, Material const & o) {
            stream << "n=" << o.real << "+i" << o.imag;
            if(o.diffuse_inclusions){

                stream << ". Diffuse inclusion: a=" << o.diffuse._albedo << ", mfp="<< o.diffuse._mean_free_path << ", fname=" << o.diffuse._fname;
            }
           
            return stream;
        }

        bool operator ==(const Material & o) const
        {
            return (real==o.real && imag==o.imag && diffuse_inclusions==o.diffuse_inclusions);
        }

        int request_n_entries(std::vector<double>& tmp){
            return (diffuse_inclusions ? diffuse.request_n_entries(tmp) : 0);
        }

        int prepare(std::vector<double>& tmp, double* storage){
            return (diffuse_inclusions ? diffuse.initialize_material(tmp, storage) : 0);
        }

        void compute_abscf(double m0real){
            double tmp = 2.0*M_PI*m0real;
            abscf = -2.0*tmp*(m0real*imag)/(m0real*m0real);
        }


    private:
        void check_arguments(double _real, double _imag){
            if(_real<=0.0) throw std::invalid_argument("non-positive refractive (real) index is not supported: ");
            if(_imag<0.0) throw std::invalid_argument("negative refractive (imag) index is not supported: ");
        }     
};

#endif