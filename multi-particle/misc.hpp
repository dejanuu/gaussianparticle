#ifndef MISC_H
#define MISC_H

#include <cmath>

inline void fill_zero(double *M, int count){
    for(int i = 0; i<count; ++i) M[i] = 0.0;
}

inline void multiply_mat_coeff(double *M, double b, int count){
    for(int i = 0; i<count; ++i) M[i] = M[i]*b;
}


inline void add_to_mat(double *M, double *A, int j, int nbins, int Asize){
    for(int i = 0; i<Asize; ++i) M[j*nbins+i] += A[i];
}

inline int convert_2d(int i, int j, int nbins){
    return i*nbins+j;
}

inline void local_to_world_space(double theta,double phi,double *X, int a, int b, int c){
    double q[3] = {X[a],X[b],X[c]};
    double mu = std::cos(theta);
    double nu = std::sqrt(1-mu*mu);
    double cphi = std::cos(phi);
    double sphi = std::sin(phi);

    X[a]=q[0]*mu*cphi-q[1]*sphi+q[2]*nu*cphi;
    X[b]=q[0]*mu*sphi+q[1]*cphi+q[2]*nu*sphi;
    X[c]=-q[0]*nu+q[2]*mu;
}



#endif