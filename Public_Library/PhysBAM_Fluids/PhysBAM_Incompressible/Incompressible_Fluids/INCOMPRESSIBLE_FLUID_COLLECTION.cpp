//#####################################################################
// Copyright 2009-2010, Nipun Kwatra, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_UNIFORM.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_FACE_ARRAYS.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Fluids/INCOMPRESSIBLE_FLUID_COLLECTION.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T_GRID> INCOMPRESSIBLE_FLUID_COLLECTION<T_GRID>::
INCOMPRESSIBLE_FLUID_COLLECTION(const T_GRID& grid_input)
    :grid(grid_input)
    //density_container(grid),temperature_container(grid)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class T_GRID> INCOMPRESSIBLE_FLUID_COLLECTION<T_GRID>::
~INCOMPRESSIBLE_FLUID_COLLECTION()
{
}
//#####################################################################
// Function Write_Output_Files
//#####################################################################
template<class T_GRID> void INCOMPRESSIBLE_FLUID_COLLECTION<T_GRID>::
Write_Output_Files(const STREAM_TYPE stream_type,const std::string& output_directory,const int frame) const
{
    std::string f=STRING_UTILITIES::string_sprintf("%d",frame);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/mac_velocities",face_velocities);
}
//#####################################################################
// Function Read_Output_Files
//#####################################################################
template<class T_GRID> void INCOMPRESSIBLE_FLUID_COLLECTION<T_GRID>::
Read_Output_Files(const STREAM_TYPE stream_type,const std::string& output_directory,const int frame)
{
    std::string f=STRING_UTILITIES::string_sprintf("%d",frame);
    std::string filename=output_directory+"/"+f+"/mac_velocities";
    std::string centered_velocity_filename=output_directory+"/"+f+"/centered_velocities";
    if(FILE_UTILITIES::File_Exists(filename)){LOG::cout<<"Reading mac_velocities "<<filename<<std::endl;
        FILE_UTILITIES::Read_From_File(stream_type,filename,face_velocities);}
    else if(FILE_UTILITIES::File_Exists(centered_velocity_filename)){
        LOG::cout<<"Converting centered velocities from "<<centered_velocity_filename<<" to mac_velocities"<<std::endl;
        T_ARRAYS_VECTOR centered_velocities;
        FILE_UTILITIES::Read_From_File(stream_type,centered_velocity_filename,centered_velocities);

        TV_INT face_index,first_cell_index,second_cell_index;int axis;
        for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){
            face_index=iterator.Face_Index();axis=iterator.Axis();
            first_cell_index=iterator.First_Cell_Index();second_cell_index=iterator.Second_Cell_Index();
            face_velocities.Component(axis)(face_index)=
                (centered_velocities(first_cell_index)(axis)+centered_velocities(second_cell_index)(axis))*(T).5;}}
}
//#####################################################################
// Function Initialize_Grids
//#####################################################################
template<class T_GRID> void INCOMPRESSIBLE_FLUID_COLLECTION<T_GRID>::
Initialize_Grids()
{
    face_velocities.Resize(grid);
}
//#####################################################################
// Function Sync_Data
//#####################################################################
template<class T_GRID> void INCOMPRESSIBLE_FLUID_COLLECTION<T_GRID>::
Sync_Data(INCOMPRESSIBLE_FLUID_COLLECTION<T_GRID>& fluid_collection,THREADED_UNIFORM_GRID<T_GRID>& threaded_grid)
{
    threaded_grid.Sync_Face_Scalar(face_velocities,fluid_collection.face_velocities);    
    threaded_grid.Sync_Scalar(viscosity,fluid_collection.viscosity);    
}
//#####################################################################
// Function Distribute_Data
//#####################################################################
template<class T_GRID> void INCOMPRESSIBLE_FLUID_COLLECTION<T_GRID>::
Distribute_Data(INCOMPRESSIBLE_FLUID_COLLECTION<T_GRID>& fluid_collection,THREADED_UNIFORM_GRID<T_GRID>& threaded_grid)
{
    threaded_grid.Distribute_Face_Scalar(face_velocities,fluid_collection.face_velocities);    
    threaded_grid.Distribute_Scalar(viscosity,fluid_collection.viscosity);    
}
//#####################################################################
template class INCOMPRESSIBLE_FLUID_COLLECTION<GRID<VECTOR<float,1> > >;
template class INCOMPRESSIBLE_FLUID_COLLECTION<GRID<VECTOR<float,2> > >;
template class INCOMPRESSIBLE_FLUID_COLLECTION<GRID<VECTOR<float,3> > >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class INCOMPRESSIBLE_FLUID_COLLECTION<GRID<VECTOR<double,1> > >;
template class INCOMPRESSIBLE_FLUID_COLLECTION<GRID<VECTOR<double,2> > >;
template class INCOMPRESSIBLE_FLUID_COLLECTION<GRID<VECTOR<double,3> > >;
#endif
