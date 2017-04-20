//#####################################################################
// Copyright 2009-2010, Nipun Kwatra, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Log/LOG.h>
#include <Grid_Tools/Grids/FACE_ITERATOR.h>
#include <Grid_Tools/Grids/GRID.h>
#include <Grid_PDE/Boundaries/BOUNDARY.h>
#include <Incompressible/Incompressible_Fluids/INCOMPRESSIBLE_FLUID_COLLECTION.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> INCOMPRESSIBLE_FLUID_COLLECTION<TV>::
INCOMPRESSIBLE_FLUID_COLLECTION(const GRID<TV>& grid_input)
    :grid(grid_input)
    //density_container(grid),temperature_container(grid)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> INCOMPRESSIBLE_FLUID_COLLECTION<TV>::
~INCOMPRESSIBLE_FLUID_COLLECTION()
{
}
//#####################################################################
// Function Write_Output_Files
//#####################################################################
template<class TV> void INCOMPRESSIBLE_FLUID_COLLECTION<TV>::
Write_Output_Files(const STREAM_TYPE stream_type,const std::string& output_directory,const int frame) const
{
    std::string f=LOG::sprintf("%d",frame);
    Write_To_File(stream_type,output_directory+"/"+f+"/mac_velocities",face_velocities);
}
//#####################################################################
// Function Read_Output_Files
//#####################################################################
template<class TV> void INCOMPRESSIBLE_FLUID_COLLECTION<TV>::
Read_Output_Files(const STREAM_TYPE stream_type,const std::string& output_directory,const int frame)
{
    std::string f=LOG::sprintf("%d",frame);
    std::string filename=output_directory+"/"+f+"/mac_velocities";
    std::string centered_velocity_filename=output_directory+"/"+f+"/centered_velocities";
    if(File_Exists(filename)){LOG::cout<<"Reading mac_velocities "<<filename<<std::endl;
        Read_From_File(stream_type,filename,face_velocities);}
    else if(File_Exists(centered_velocity_filename)){
        LOG::cout<<"Converting centered velocities from "<<centered_velocity_filename<<" to mac_velocities"<<std::endl;
        ARRAY<TV,TV_INT> centered_velocities;
        Read_From_File(stream_type,centered_velocity_filename,centered_velocities);

        TV_INT face_index,first_cell_index,second_cell_index;int axis;
        for(FACE_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()){
            face_index=iterator.Face_Index();axis=iterator.Axis();
            first_cell_index=iterator.First_Cell_Index();second_cell_index=iterator.Second_Cell_Index();
            face_velocities.Component(axis)(face_index)=
                (centered_velocities(first_cell_index)(axis)+centered_velocities(second_cell_index)(axis))*(T).5;}}
}
//#####################################################################
// Function Initialize_Grids
//#####################################################################
template<class TV> void INCOMPRESSIBLE_FLUID_COLLECTION<TV>::
Initialize_Grids()
{
    face_velocities.Resize(grid);
}
//#####################################################################
namespace PhysBAM{
template class INCOMPRESSIBLE_FLUID_COLLECTION<VECTOR<float,1> >;
template class INCOMPRESSIBLE_FLUID_COLLECTION<VECTOR<float,2> >;
template class INCOMPRESSIBLE_FLUID_COLLECTION<VECTOR<float,3> >;
template class INCOMPRESSIBLE_FLUID_COLLECTION<VECTOR<double,1> >;
template class INCOMPRESSIBLE_FLUID_COLLECTION<VECTOR<double,2> >;
template class INCOMPRESSIBLE_FLUID_COLLECTION<VECTOR<double,3> >;
}
