//#####################################################################
// Copyright 2009, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Arrays_Nd/ARRAYS_ND.h>
#include <Core/Utilities/VIEWER_DIR.h>
#include <Grid_Tools/Grids/GRID.h>
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
Write_Output_Files(const STREAM_TYPE stream_type,const VIEWER_DIR& viewer_dir) const
{
    Write_To_File(stream_type,viewer_dir.current_directory+"/psi",psi);
    Write_To_File(stream_type,viewer_dir.current_directory+"/euler_U",U);

    // TODO(jontg): Write this out optionally.
    COMPRESSIBLE_AUXILIARY_DATA::Write_Auxiliary_Files(stream_type,viewer_dir,*this,false);
}
//#####################################################################
// Function Read_Output_Files
//#####################################################################
template<class TV> void COMPRESSIBLE_FLUID_COLLECTION<TV>::
Read_Output_Files(const VIEWER_DIR& viewer_dir)
{
    if(File_Exists(viewer_dir.current_directory+"/psi")){
        Read_From_File(viewer_dir.current_directory+"/psi",psi);}

    if(File_Exists(viewer_dir.current_directory+"/euler_U")){
        Read_From_File(viewer_dir.current_directory+"/euler_U",U);}
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
