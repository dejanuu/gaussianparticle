#ifndef MATERIALS_H
#define MATERIALS_H


#include "inputreader.hpp"
#include "definitions.hpp"
#include "inputreader.hpp"
#include "smaterial.hpp"


#include <sstream>
#include <iostream>
#include <stdexcept>
#include <fstream>
#include <map>
#include <tuple>
#include <list>

class Materials{
    public:
        Materials(InputReader& reader) : materials(){
            bool contains_diffuse_inclusions = false;
            std::cout << create_heading("MATERIALS",HEADINGS::HEADING) << std::endl;
            std::cout << "Material format: ref_real ref_imag diffuse_inclusions mean_free_path albedo fname"  << std::endl;
            try{
                std::cout << "Surrouding medium format: ref_real ref_imag 0" << std::endl;
                materials.push_back(reader.extract<Material>("media",Material(1.0,0.0)));
                materials.back().compute_abscf(materials[0].real);
                if(materials[0].diffuse_inclusions) throw std::logic_error("Surrouding medium doesn't support diffuse inclusions");
                for(int i=1;i<255;++i){
                    std::string material = std::string("material") + std::to_string(i);
                    if(reader.arg_found(material)){
                        materials.push_back(reader.extract<Material>(material,Material(1.3,0.0)));
                        materials.back().compute_abscf(materials[0].real); 
                        if(materials.back().diffuse_inclusions) contains_diffuse_inclusions = true; 
                    }
                }
                
                materials_tmp.insert(std::pair<int,Material*>(0,&materials[0]));
                int indx=0;
                for(int i=1;i<255;++i){
                    std::string material = std::string("material") + std::to_string(i);
                    if(reader.arg_found(material)){
                        materials_tmp.insert(std::pair<int,Material*>(i,&materials[++indx]));
                    }
                }
                if(contains_diffuse_inclusions){
                    std::cout << create_heading("CREATING DIFFUSE SCATTERERS",HEADINGS::SUB_HEADING) << std::endl;
                    initialize_materials();
                }
                success = true;
            }catch (const std::exception& e) {
                std::cout << create_heading("#!#!#!#ERROR IN MATERIALS#!#!#!#",HEADINGS::HEADING) << std::endl;
                std::cout << e.what() << std::endl;
                return;
            }
            
            std::cout << create_heading("",HEADINGS::HEADING) << std::endl;
        }

        bool connect_material(std::vector<std::tuple<int, int>>& material_inds, std::list<Triangle>& triangles){
            if(!success) return false;
            try{
                int i = 0;
                for(auto it = triangles.begin(); it!=triangles.end();++it){
                    //!Add error checking
                    //if(materials_initial.find(std::get<0>(*it))!=materials_initial.end()){
                    //    materials.
                    //}else{
                    //    throw std::invalid_argument(std::string("Material not given: ")+std::to_string(std::get<0>(*it)));
                    //}

                    Material * mat1=materials_tmp[std::get<0>(material_inds[i])];
                    Material * mat2=materials_tmp[std::get<1>(material_inds[i])];
                    materials_storage.insert(std::pair<Triangle*,std::tuple<Material*,Material*>>
                                    (&*it,std::tuple<Material*, Material*>(mat1,mat2)));
                    ++i;
                }
            }catch (const std::exception& e) {
                std::cout << e.what() << std::endl;
                return false;
            }
            return true;
        }

        void initialize_materials(){
            int size = 0;
            for(auto material= materials.begin();material!=materials.end();++material){
                size += material->request_n_entries(tmp_storage);
            }
            storage = new double[size];
            int indx=0;
            for(auto material= materials.begin();material!=materials.end();++material){
                
                indx += material->prepare(tmp_storage,&storage[indx]);
            }
        }

        ~Materials(){
            delete[] storage;
        }




        Material* get_first_material(Triangle* triangle){
            return std::get<0>(materials_storage[triangle]);
        }

        Material* get_second_material(Triangle* triangle){
            return std::get<1>(materials_storage[triangle]);
        }

        Material* get_surrounding_material(){
            return &materials[0];
        }

    private:
        double* storage = nullptr;
        int entries;
        std::vector<Material> materials;
        std::map<int,Material*> materials_tmp;
        std::map<Triangle*,std::tuple<Material*,Material*>> materials_storage;
        bool success = false;
        std::string empty = "null";
        std::vector<double> tmp_storage;
};


#endif
