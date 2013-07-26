//#####################################################################
// Copyright 2009, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Grids_Uniform/GRID.h>
#include <Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <Compressible/Compressible_Fluids/COMPRESSIBLE_AUXILIARY_DATA.h>
#include <Compressible/Compressible_Fluids/COMPRESSIBLE_FLUID_COLLECTION.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> COMPRESSIBLE_FLUID_COLLECTION<TV>::
COMPRESSIBLE_FLUID_COLLECTION(const GRID<TV>& grid_input)
:grid(grid_input),eos(0),U()
{}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> COMPRESSIBLE_FLUID_COLLECTION<TV>::
~COMPRESSIBLE_FLUID_COLLECTION()
{}
//#####################################################################
// Function Write_Output_Files
//#####################################################################
template<class TV> void COMPRESSIBLE_FLUID_COLLECTION<TV>::
Write_Output_Files(const STREAM_TYPE stream_type,const std::string& output_directory,const int frame) const
{
    std::string f=STRING_UTILITIES::string_sprintf("%d",frame);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/psi",psi);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/euler_U",U);

    // TODO(jontg): Write this out optionally.
    COMPRESSIBLE_AUXILIARY_DATA::Write_Auxiliary_Files(stream_type,output_directory+"/"+f,frame,*this,false);
}
//#####################################################################
// Function Read_Output_Files
//#####################################################################
template<class TV> void COMPRESSIBLE_FLUID_COLLECTION<TV>::
Read_Output_Files(const STREAM_TYPE stream_type,const std::string& output_directory,const int frame)
{
    std::string f=STRING_UTILITIES::string_sprintf("%d",frame);
    if(FILE_UTILITIES::File_Exists(output_directory+"/"+f+"/psi")){
        FILE_UTILITIES::Read_From_File(stream_type,output_directory+"/"+f+"/psi",psi);}

    if(FILE_UTILITIES::File_Exists(output_directory+"/"+f+"/euler_U")){
        FILE_UTILITIES::Read_From_File(stream_type,output_directory+"/"+f+"/euler_U",U);}
}
//#####################################################################
// Function Initialize_Grids
//#####################################################################
template<class TV> void COMPRESSIBLE_FLUID_COLLECTION<TV>::
Initialize_Grids()
{
    U.Resize(grid.Domain_Indices());
}
//#####################################################################
template class COMPRESSIBLE_FLUID_COLLECTION<VECTOR<float,1> >;
template class COMPRESSIBLE_FLUID_COLLECTION<VECTOR<float,2> >;
template class COMPRESSIBLE_FLUID_COLLECTION<VECTOR<float,3> >;
template class COMPRESSIBLE_FLUID_COLLECTION<VECTOR<double,1> >;
template class COMPRESSIBLE_FLUID_COLLECTION<VECTOR<double,2> >;
template class COMPRESSIBLE_FLUID_COLLECTION<VECTOR<double,3> >;
}
