//#####################################################################
// Copyright 2009, Michael Lentine, Avi Robinson-Mosher, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Utilities/PROCESS_UTILITIES.h>
#include <PhysBAM_Geometry/Geometry_Particles/REGISTER_GEOMETRY_READ_WRITE.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_RIGID_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Forces/FLUID_GRAVITY.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Forces/INCOMPRESSIBILITY.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/INCOMPRESSIBLE_EXAMPLE.h>
#include <PhysBAM_Dynamics/Advection_Equations/ADVECTION_CONSERVATIVE_UNIFORM.h>
#include <PhysBAM_Dynamics/Advection_Equations/ADVECTION_CONSERVATIVE_UNIFORM_FORWARD.h>
using namespace PhysBAM;
//#####################################################################
// INCOMPRESSIBLE_EXAMPLE
//#####################################################################
template<class TV> INCOMPRESSIBLE_EXAMPLE<TV>::
INCOMPRESSIBLE_EXAMPLE(const STREAM_TYPE stream_type_input)
    :stream_type(stream_type_input),initial_time(0),first_frame(0),last_frame(100),frame_rate(24),
    restart(0),write_debug_data(false),analytic_test(false),order(1),output_directory("output"),
    number_of_ghost_cells(3),cfl((T).9),use_viscosity(false),mac_grid(TV_INT(),RANGE<TV>::Unit_Box(),true),mpi_grid(0),//incompressible_fluid_collection(mac_grid),
    projection(mac_grid),incompressible(mac_grid,projection),boundary(0),rigid_geometry_collection(this)
{
    PROCESS_UTILITIES::Set_Floating_Point_Exception_Handling(true);
    PROCESS_UTILITIES::Set_Backtrace(true);
    incompressible.Set_Custom_Advection(advection_scalar);
    //incompressible.Add_Force(new INCOMPRESSIBILITY<GRID<TV> >(projection));
    //incompressible.Add_Force(new FLUID_GRAVITY<GRID<TV> >());
    for(int i=0;i<TV::dimension;i++){domain_boundary(i)(0)=true;domain_boundary(i)(1)=true;}
    Initialize_Geometry_Particle();
}
//#####################################################################
// ~INCOMPRESSIBLE_EXAMPLE
//#####################################################################
template<class TV> INCOMPRESSIBLE_EXAMPLE<TV>::
~INCOMPRESSIBLE_EXAMPLE()
{
    if(mpi_grid) delete boundary;
}
//#####################################################################
// 
//#####################################################################
template<class TV> void INCOMPRESSIBLE_EXAMPLE<TV>::
Write_Output_Files(const int frame)
{
    std::string f=STRING_UTILITIES::string_sprintf("%d",frame);
    if(mpi_grid) FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/common/global_grid",mpi_grid->global_grid);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/grid",mac_grid);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/common/grid",mac_grid);
    if(incompressible.use_analytic_energy) FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/analytic_energy",incompressible.analytic_energy);
    ADVECTION_CONSERVATIVE_UNIFORM<GRID<TV>,T>* advection_conservative=dynamic_cast<ADVECTION_CONSERVATIVE_UNIFORM<GRID<TV>,T>*>(incompressible.advection);
    if(advection_conservative){
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/weights_barjc",advection_conservative->sum_jc);
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/weights_barjc_cell",advection_conservative->sum_jc_cell);}
    if(write_debug_data){
        if(advection_conservative){
            FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/weights_i",advection_conservative->weights_to);
            FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/weights_j",advection_conservative->weights_from);}
        ARRAY<T,FACE_INDEX<TV::dimension> > face_velocities_ghost(mac_grid,number_of_ghost_cells,false);
        ARRAY<T,TV_INT> density_ghost(mac_grid.Domain_Indices(number_of_ghost_cells),false);
        incompressible.boundary->Fill_Ghost_Cells_Face(mac_grid,face_velocities,face_velocities_ghost,0,number_of_ghost_cells);
        incompressible.boundary->Fill_Ghost_Cells(mac_grid,density,density_ghost,(T)0,0,number_of_ghost_cells);
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/mac_velocities",face_velocities_ghost);
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/density",density_ghost);
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/pressure",incompressible.projection.p);
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/psi_N",projection.elliptic_solver->psi_N);
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/psi_D",projection.elliptic_solver->psi_D);}
    else{
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/mac_velocities",face_velocities);
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/density",density);}
    if(!stream_type.use_doubles)
        Read_Write<RIGID_GEOMETRY_COLLECTION<TV>,float>::Write(stream_type,output_directory,frame,rigid_geometry_collection);
    else
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
        Read_Write<RIGID_GEOMETRY_COLLECTION<TV>,double>::Write(stream_type,output_directory,frame,rigid_geometry_collection);
#else
        PHYSBAM_FATAL_ERROR("Cannot read doubles");
#endif
}
template<class TV> void INCOMPRESSIBLE_EXAMPLE<TV>::
Read_Output_Files(const int frame)
{
    std::string f=STRING_UTILITIES::string_sprintf("%d",frame);
    FILE_UTILITIES::Read_From_File(stream_type,output_directory+"/"+f+"/density",density);
    if(incompressible.use_analytic_energy) FILE_UTILITIES::Read_From_File(stream_type,output_directory+"/"+f+"/analytic_energy",incompressible.analytic_energy);
    std::string filename;
    filename=output_directory+"/"+f+"/mac_velocities";
    if(FILE_UTILITIES::File_Exists(filename)){LOG::cout<<"Reading mac_velocities "<<filename<<std::endl;FILE_UTILITIES::Read_From_File(stream_type,filename,face_velocities);}
    filename=output_directory+"/"+f+"/pressure";
    if(FILE_UTILITIES::File_Exists(filename)){LOG::cout<<"Reading pressure "<<filename<<std::endl;FILE_UTILITIES::Read_From_File(stream_type,filename,incompressible.projection.p);}
    ADVECTION_CONSERVATIVE_UNIFORM<GRID<TV>,T>* advection_conservative=dynamic_cast<ADVECTION_CONSERVATIVE_UNIFORM<GRID<TV>,T>*>(incompressible.advection);
    if(advection_conservative){
        filename=output_directory+"/"+f+"/weights_barjc";
        if(FILE_UTILITIES::File_Exists(filename)){LOG::cout<<"Reading wjc "<<filename<<std::endl;FILE_UTILITIES::Read_From_File(stream_type,filename,advection_conservative->sum_jc);}
        filename=output_directory+"/"+f+"/weights_barjc_cell";
        if(FILE_UTILITIES::File_Exists(filename)){LOG::cout<<"Reading wjc cell "<<filename<<std::endl;FILE_UTILITIES::Read_From_File(stream_type,filename,advection_conservative->sum_jc_cell);}}
    if(!stream_type.use_doubles)
        Read_Write<RIGID_GEOMETRY_COLLECTION<TV>,float>::Read(stream_type,output_directory,frame,rigid_geometry_collection);
    else
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
        Read_Write<RIGID_GEOMETRY_COLLECTION<TV>,double>::Read(stream_type,output_directory,frame,rigid_geometry_collection);
#else
        PHYSBAM_FATAL_ERROR("Cannot read doubles");
#endif
}
//#####################################################################
template class INCOMPRESSIBLE_EXAMPLE<VECTOR<float,1> >;
template class INCOMPRESSIBLE_EXAMPLE<VECTOR<float,2> >;
template class INCOMPRESSIBLE_EXAMPLE<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class INCOMPRESSIBLE_EXAMPLE<VECTOR<double,1> >;
template class INCOMPRESSIBLE_EXAMPLE<VECTOR<double,2> >;
template class INCOMPRESSIBLE_EXAMPLE<VECTOR<double,3> >;
#endif
