//#####################################################################
// Copyright 2009-2010, Michael Lentine, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
//#include <PhysBAM_Tools/Grids_Uniform_Advection/ADVECTION_SEMI_LAGRANGIAN_UNIFORM_BETA.h>
#include <PhysBAM_Tools/Parallel_Computation/DOMAIN_ITERATOR_THREADED.h>
#include "ADVECTION_EXAMPLE.h"
#include <pthread.h>
using namespace PhysBAM;
//#####################################################################
// ADVECTION_EXAMPLE
//#####################################################################
template<class TV> ADVECTION_EXAMPLE<TV>::
ADVECTION_EXAMPLE(const STREAM_TYPE stream_type_input,int number_of_threads)
    :stream_type(stream_type_input),initial_time(0),first_frame(0),last_frame(100),frame_rate(24),
    restart(0),write_debug_data(false),output_directory("output"),cfl(.9),thread_queue(number_of_threads>1?new THREAD_QUEUE(number_of_threads):0),
    mac_grid(TV_INT(),RANGE<TV>::Unit_Box(),true),mpi_grid(0),advection_scalar(thread_queue),boundary(0)
{
    pthread_mutex_init(&lock,0);
}
//#####################################################################
// ~ADVECTION_EXAMPLE
//#####################################################################
template<class TV> ADVECTION_EXAMPLE<TV>::
~ADVECTION_EXAMPLE()
{
    if(mpi_grid) delete boundary;
}
//#####################################################################
// CFL 
//#####################################################################
template<class TV> typename TV::SCALAR ADVECTION_EXAMPLE<TV>::
CFL(T_FACE_ARRAYS_SCALAR& face_velocities)
{
    T dt=FLT_MAX;
    DOMAIN_ITERATOR_THREADED_ALPHA<ADVECTION_EXAMPLE<TV>,TV>(mac_grid.Domain_Indices(),thread_queue).template Run<T_FACE_ARRAYS_SCALAR&,T&>(*this,&ADVECTION_EXAMPLE::CFL_Threaded,face_velocities,dt);
    return dt;
}
//#####################################################################
// CFL_Threaded 
//#####################################################################
template<class TV> void ADVECTION_EXAMPLE<TV>::
CFL_Threaded(RANGE<TV_INT>& domain,T_FACE_ARRAYS_SCALAR& face_velocities,T& dt)
{
    T dt_convection=0;
    for(CELL_ITERATOR iterator(mac_grid,domain);iterator.Valid();iterator.Next()){
        TV_INT cell=iterator.Cell_Index();T local_V_norm=0;
        for(int axis=0;axis<GRID<TV>::dimension;axis++)
            local_V_norm+=mac_grid.one_over_dX[axis]*maxabs(face_velocities(axis,mac_grid.First_Face_Index_In_Cell(axis,cell)),face_velocities(axis,mac_grid.Second_Face_Index_In_Cell(axis,cell)));
        dt_convection=max(dt_convection,local_V_norm);}
    pthread_mutex_lock(&lock);
    dt=min(dt,(T)1.0/dt_convection);
    pthread_mutex_unlock(&lock);
}
//#####################################################################
// Write_Output_Files 
//#####################################################################
template<class TV> void ADVECTION_EXAMPLE<TV>::
Write_Output_Files(const int frame)
{
    std::string f=STRING_UTILITIES::string_sprintf("%d",frame);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/mac_velocities",face_velocities);
    if(mpi_grid) FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/common/global_grid",mpi_grid->global_grid);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/grid",mac_grid);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/common/grid",mac_grid);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/density",density);
}
template<class TV> void ADVECTION_EXAMPLE<TV>::
Read_Output_Files(const int frame)
{
    std::string f=STRING_UTILITIES::string_sprintf("%d",frame);
    FILE_UTILITIES::Read_From_File(stream_type,output_directory+"/"+f+"/density",density);
    std::string filename;
    filename=output_directory+"/"+f+"/mac_velocities";
    if(FILE_UTILITIES::File_Exists(filename)){LOG::cout<<"Reading mac_velocities "<<filename<<std::endl;FILE_UTILITIES::Read_From_File(stream_type,filename,face_velocities);}
}
//#####################################################################
template class ADVECTION_EXAMPLE<VECTOR<float,1> >;
template class ADVECTION_EXAMPLE<VECTOR<float,2> >;
template class ADVECTION_EXAMPLE<VECTOR<float,3> >;
template class ADVECTION_EXAMPLE<VECTOR<double,1> >;
template class ADVECTION_EXAMPLE<VECTOR<double,2> >;
template class ADVECTION_EXAMPLE<VECTOR<double,3> >;
