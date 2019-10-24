# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( SIRISn CXX C Fortran)

cmake_minimum_required(VERSION 2.8.10)

find_package(CGAL REQUIRED)
find_package(Boost REQUIRED)


  

if ( CGAL_FOUND )
    include( ${CGAL_USE_FILE} )


    set(CMAKE_Fortran_FLAGS "-O3 -cpp -mtune=native -std=f2008 -ffree-form")
    set(CMAKE_CXX_FLAGS "-O3 -std=c++11 -mtune=native -frounding-math")
    #set(CMAKE_CXX_FLAGS_RELEASE "-g -frounding-math -std=c++11")



    #create_single_source_cgal_program(  )
    add_executable(SIRISn "inputreader.cpp" "misc.hpp" "rng.cpp" 
                        "offreader.cpp" "meshreader.cpp" "materials.cpp"
                        "mray.cpp" "sray.cpp" "detector.cpp" "physicsengine.cpp" 
                        "smaterial.cpp" "geometry.cpp" "outputwriter.cpp"
                        "physicsengine.cpp" "tracer.cpp" "definitions.hpp" 
                        "sirisinterface.f90" "../src/sirisconstants.f" 
                        "splinetools.f90" "mathroutines.f90" "../src/sirisradtrans.f"
                        "../src/sirismath.f" "../src/sirismaterial.f"
                        "../src/sirisnumint.f" "../src/sirisray.f"
                        "main.cpp")
    target_link_libraries(SIRISn ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})


else()

    message(STATUS "CGAL missing!!!")

endif()

