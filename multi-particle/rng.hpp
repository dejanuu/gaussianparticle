#ifndef RNG_H
#define RNG_H

#include "definitions.hpp"
#include <random>
#include <algorithm>




class RNG{
    public:
        //valgrind doesn't like this
        //RNG(int seed) : _seed{seed<=0 ? 1 : _seed}, dist{0.0,1.0}, twister{_seed}{
        RNG(int seed) : _seed{seed<=0 ? (int)std::random_device()() : seed}, dist{0.0,1.0}, twister{_seed}{
            std::cout << "Using seed: " << _seed << std::endl;
        };
    
        inline double rand(){
            return dist(twister);
        }

        inline void rand(std::vector<double>& vec, int num_items){
            auto rng = std::bind(std::ref(dist), std::ref(twister));
            std::generate_n(vec.begin(), num_items, rng);
            //std::cout << vec[0] << std::endl;
        }

    private:
        int _seed;
        std::uniform_real_distribution<double> dist;
        std::mt19937_64 twister;
        
};


#endif