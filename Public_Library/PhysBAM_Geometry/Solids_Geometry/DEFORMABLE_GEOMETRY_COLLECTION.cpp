//#####################################################################
// Copyright 2003-2009, Zhaosheng Bao, Ron Fedkiw, Geoffrey Irving, Sergey Levine, Igor Neverov, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DEFORMABLE_GEOMETRY_COLLECTION
//#####################################################################
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Read_Write/FILE_UTILITIES.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Solids_Geometry/DEFORMABLE_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/STRUCTURE.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> DEFORMABLE_GEOMETRY_COLLECTION<TV>::
DEFORMABLE_GEOMETRY_COLLECTION(GEOMETRY_PARTICLES<TV>& particles_input)
    :particles(particles_input)
{
    particles.Store_Velocity();
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> DEFORMABLE_GEOMETRY_COLLECTION<TV>::
~DEFORMABLE_GEOMETRY_COLLECTION()
{
    structures.Delete_Pointers_And_Clean_Memory();
}
//#####################################################################
// Function Read
//#####################################################################
template<class TV> void DEFORMABLE_GEOMETRY_COLLECTION<TV>::
Read(const STREAM_TYPE stream_type,const std::string& prefix,const std::string& static_prefix,const int frame,const int static_frame,const bool include_static_variables)
{
    Read_Dynamic_Variables(stream_type,prefix,frame);
    if(include_static_variables) Read_Static_Variables(stream_type,static_prefix,static_frame);
}
//#####################################################################
// Function Write
//#####################################################################
template<class TV> void DEFORMABLE_GEOMETRY_COLLECTION<TV>::
Write(const STREAM_TYPE stream_type,const std::string& prefix,const std::string& static_prefix,const int frame,const int static_frame,const bool include_static_variables) const
{
    Write_Dynamic_Variables(stream_type,prefix,frame);
    if(include_static_variables) Write_Static_Variables(stream_type,static_prefix,static_frame);
}
//#####################################################################
// Function Read_Static_Variables
//#####################################################################
template<class TV> void DEFORMABLE_GEOMETRY_COLLECTION<TV>::
Read_Static_Variables(const STREAM_TYPE stream_type,const std::string& prefix,const int frame)
{
    std::string f=frame==-1?"common":FILE_UTILITIES::Number_To_String(frame);
    std::istream* input_raw=FILE_UTILITIES::Safe_Open_Input(prefix+"/"+f+"/deformable_object_structures");
    TYPED_ISTREAM input(*input_raw,stream_type);
    int m;Read_Binary(input,m);
    // TODO: merge this functionality with dynamic lists to allow for more flexibility
    if(!structures.m){ // // create and read all structures from scratch
        structures.Resize(m);
        for(int k=0;k<structures.m;k++) structures(k)=STRUCTURE<TV>::Create_Structure(input,particles);
    }
    else if(structures.m<=m){
        int old_number_of_structures=structures.m;structures.Resize(m);
        for(int k=0;k<old_number_of_structures;k++) structures(k)->Read_Structure(input);
        for(int k=old_number_of_structures;k<m;k++) structures(k)=STRUCTURE<TV>::Create_Structure(input,particles);}
    else{
        LOG::cout<<"Current number of structures ("<<structures.m<<") is greater than number in file ("<<m<<").";
        PHYSBAM_FATAL_ERROR();}
    delete input_raw;
}
//#####################################################################
// Function Read_Dynamic_Variables
//#####################################################################
template<class TV> void DEFORMABLE_GEOMETRY_COLLECTION<TV>::
Read_Dynamic_Variables(const STREAM_TYPE stream_type,const std::string& prefix,const int frame)
{
    FILE_UTILITIES::Read_From_File(stream_type,prefix+"/"+FILE_UTILITIES::Number_To_String(frame)+"/deformable_object_particles",particles);
    // if number==0, the particles format doesn't remember the set of attributes, so the following line makes restarts look more exact
    particles.Store_Velocity();
}
//#####################################################################
// Function Write_Static_Variables
//#####################################################################
template<class TV> void DEFORMABLE_GEOMETRY_COLLECTION<TV>::
Write_Static_Variables(const STREAM_TYPE stream_type,const std::string& prefix,const int frame) const
{
    std::string f=frame==-1?"common":FILE_UTILITIES::Number_To_String(frame);
    std::ostream* output_raw=FILE_UTILITIES::Safe_Open_Output(prefix+"/"+f+"/deformable_object_structures");
    TYPED_OSTREAM output(*output_raw,stream_type);
    Write_Binary(output,structures.m);
    for(int k=0;k<structures.m;k++) structures(k)->Write_Structure(output);
    delete output_raw;
}
//#####################################################################
// Function Write_Dynamic_Variables
//#####################################################################
template<class TV> void DEFORMABLE_GEOMETRY_COLLECTION<TV>::
Write_Dynamic_Variables(const STREAM_TYPE stream_type,const std::string& prefix,const int frame) const
{
    FILE_UTILITIES::Write_To_File(stream_type,prefix+"/"+FILE_UTILITIES::Number_To_String(frame)+"/deformable_object_particles",particles);
}
//#####################################################################
template class DEFORMABLE_GEOMETRY_COLLECTION<VECTOR<float,1> >;
template class DEFORMABLE_GEOMETRY_COLLECTION<VECTOR<float,2> >;
template class DEFORMABLE_GEOMETRY_COLLECTION<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class DEFORMABLE_GEOMETRY_COLLECTION<VECTOR<double,1> >;
template class DEFORMABLE_GEOMETRY_COLLECTION<VECTOR<double,2> >;
template class DEFORMABLE_GEOMETRY_COLLECTION<VECTOR<double,3> >;
#endif

