//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Grids_Uniform/GRID.h>
#include <Incompressible/Incompressible_Fluids/INCOMPRESSIBLE_FLUID_COLLECTION.h>
#include <Fluids/Fluids/FLUID_COLLECTION.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> FLUID_COLLECTION<TV>::
FLUID_COLLECTION(const GRID<TV>& grid_input)
    :grid(grid_input),incompressible_fluid_collection(grid)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> FLUID_COLLECTION<TV>::
~FLUID_COLLECTION()
{
}
//#####################################################################
// Function Write_Output_Files
//#####################################################################
template<class TV> void FLUID_COLLECTION<TV>::
Write_Output_Files(const STREAM_TYPE stream_type,const std::string& output_directory,const int frame) const
{
    //compressible_fluid_collection.Write_Output_Files(frame);
    incompressible_fluid_collection.Write_Output_Files(stream_type,output_directory,frame);
}
//#####################################################################
// Function Read_Output_Files
//#####################################################################
template<class TV> void FLUID_COLLECTION<TV>::
Read_Output_Files(const STREAM_TYPE stream_type,const std::string& output_directory,const int frame)
{
    //compressible_fluid_collection.Read_Output_Files(frame);
    incompressible_fluid_collection.Read_Output_Files(stream_type,output_directory,frame);
}
//#####################################################################
// Function Initialize
//#####################################################################
template<class TV> void FLUID_COLLECTION<TV>::
Initialize_Grids()
{
    //compressible_fluid_collection.Initialize_Grids();
    incompressible_fluid_collection.Initialize_Grids();
}
//#####################################################################
namespace PhysBAM{
template class FLUID_COLLECTION<VECTOR<float,1> >;
template class FLUID_COLLECTION<VECTOR<float,2> >;
template class FLUID_COLLECTION<VECTOR<float,3> >;
template class FLUID_COLLECTION<VECTOR<double,1> >;
template class FLUID_COLLECTION<VECTOR<double,2> >;
template class FLUID_COLLECTION<VECTOR<double,3> >;
}

