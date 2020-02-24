#ifndef OFFREADER_H
#define OFFREADER_H

#include "definitions.hpp"
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <fstream>
#include <string>
#include <vector>
#include <tuple>
#include <boost/algorithm/string.hpp>


class OFFReader{
    public:
        void read_file(std::string& fname, std::vector<Point>& vertices, 
                                std::list<Triangle>& triangles, 
                                std::vector<std::tuple<int, int>>& material_inds, 
                                bool centralize_vertices);
};


#endif
