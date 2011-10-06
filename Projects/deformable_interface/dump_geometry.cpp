//#####################################################################
// Copyright 2007, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Geometry/Basic_Geometry/BASIC_GEOMETRY_FORWARD.h>
#include <PhysBAM_Geometry/Geometry_Particles/REGISTER_GEOMETRY_READ_WRITE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <fstream>
#include "../data_exchange/deformable_body_simulation.h"
#include "DATA_EXCHANGE_CONVERSION.h"
#include "DEFORMABLE_EXAMPLE.h"
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
using namespace PhysBAM;

int main(int argc,char* argv[])
{
    typedef float T;
    typedef float RW;
    STREAM_TYPE stream_type((RW()));
    typedef VECTOR<T,3> TV;

    Initialize_Geometry_Particle();
    Initialize_Read_Write_Structures();

    std::string scene_file=argv[1];
    int read_text=atoi(argv[2]);

    data_exchange::deformable_body_simulation dbs;
    std::ifstream ifs(scene_file.c_str());
    if(read_text){
        boost::archive::text_iarchive ia(ifs);
        ia >> dbs;}
    else{
        boost::archive::binary_iarchive ia(ifs);
        ia >> dbs;}

    for(size_t i=0; i<dbs.simulation_objects.size(); i++){
        TRIANGLE_MESH tm;
        GEOMETRY_PARTICLES<TV> particles;
        TRIANGULATED_SURFACE<T> surface(tm,particles);
        if(data_exchange::deformable_body* body=dynamic_cast<data_exchange::deformable_body*>(dbs.simulation_objects[i])){
            Triangulated_Surface_From_Data_Exchange(surface,body->mesh,body->position,0);}
        else if(data_exchange::scripted_geometry* body=dynamic_cast<data_exchange::scripted_geometry*>(dbs.simulation_objects[i])){
            Triangulated_Surface_From_Data_Exchange(surface,body->mesh,body->position,0);}
        else continue;
        FILE_UTILITIES::Write_To_File(stream_type,STRING_UTILITIES::string_sprintf("%s.geo.%d.tri",scene_file.c_str(),i),surface);}

    return 0;
}
