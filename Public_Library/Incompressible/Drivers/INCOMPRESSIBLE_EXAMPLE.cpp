//#####################################################################
// Copyright 2009, Michael Lentine, Avi Robinson-Mosher, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Utilities/PROCESS_UTILITIES.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <Incompressible/Drivers/INCOMPRESSIBLE_EXAMPLE.h>
#include <Incompressible/Forces/FLUID_GRAVITY.h>
#include <Incompressible/Forces/INCOMPRESSIBILITY.h>
using namespace PhysBAM;
//#####################################################################
// INCOMPRESSIBLE_EXAMPLE
//#####################################################################
template<class TV> INCOMPRESSIBLE_EXAMPLE<TV>::
INCOMPRESSIBLE_EXAMPLE(const STREAM_TYPE stream_type_input)
    :stream_type(stream_type_input),last_frame(100),frame_rate(24),
    restart(0),write_debug_data(false),analytic_test(false),order(1),viewer_dir("output"),
    number_of_ghost_cells(3),cfl((T).9),use_viscosity(false),mac_grid(TV_INT(),RANGE<TV>::Unit_Box(),true),mpi_grid(0),//incompressible_fluid_collection(mac_grid),
    projection(mac_grid),incompressible(mac_grid,projection),boundary(0),rigid_body_collection(0)
{
    PROCESS_UTILITIES::Set_Floating_Point_Exception_Handling(true);
    PROCESS_UTILITIES::Set_Backtrace(true);
    incompressible.Set_Custom_Advection(advection_scalar);
    //incompressible.Add_Force(new INCOMPRESSIBILITY<TV>(projection));
    //incompressible.Add_Force(new FLUID_GRAVITY<TV>());
    for(int i=0;i<TV::m;i++){domain_boundary(i)(0)=true;domain_boundary(i)(1)=true;}
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
Write_Output_Files()
{
    if(mpi_grid) Write_To_File(stream_type,viewer_dir.output_directory+"/common/global_grid",mpi_grid->global_grid);
    Write_To_File(stream_type,viewer_dir.current_directory+"/grid",mac_grid);
    Write_To_File(stream_type,viewer_dir.output_directory+"/common/grid",mac_grid);
    if(write_debug_data){
        ARRAY<T,FACE_INDEX<TV::m> > face_velocities_ghost(mac_grid,number_of_ghost_cells,no_init);
        ARRAY<T,TV_INT> density_ghost(mac_grid.Domain_Indices(number_of_ghost_cells),no_init);
        incompressible.boundary->Fill_Ghost_Faces(mac_grid,face_velocities,face_velocities_ghost,0,number_of_ghost_cells);
        incompressible.boundary->Fill_Ghost_Cells(mac_grid,density,density_ghost,(T)0,0,number_of_ghost_cells);
        Write_To_File(stream_type,viewer_dir.current_directory+"/mac_velocities",face_velocities_ghost);
        Write_To_File(stream_type,viewer_dir.current_directory+"/density",density_ghost);
        Write_To_File(stream_type,viewer_dir.current_directory+"/pressure",incompressible.projection.p);
        Write_To_File(stream_type,viewer_dir.current_directory+"/psi_N",projection.elliptic_solver->psi_N);
        Write_To_File(stream_type,viewer_dir.current_directory+"/psi_D",projection.elliptic_solver->psi_D);}
    else{
        Write_To_File(stream_type,viewer_dir.current_directory+"/mac_velocities",face_velocities);
        Write_To_File(stream_type,viewer_dir.current_directory+"/density",density);}
    rigid_body_collection.Write(stream_type,viewer_dir);
    PHYSBAM_FATAL_ERROR("Cannot read doubles");
}
template<class TV> void INCOMPRESSIBLE_EXAMPLE<TV>::
Read_Output_Files()
{
    Read_From_File(viewer_dir.current_directory+"/density",density);
    std::string filename;
    filename=viewer_dir.current_directory+"/mac_velocities";
    if(File_Exists(filename)){
        LOG::cout<<"Reading mac_velocities "<<filename<<std::endl;
        Read_From_File(filename,face_velocities);}
    filename=viewer_dir.current_directory+"/pressure";
    if(File_Exists(filename)){
        LOG::cout<<"Reading pressure "<<filename<<std::endl;
        Read_From_File(filename,incompressible.projection.p);}
    rigid_body_collection.Read(viewer_dir);
}
//#####################################################################
namespace PhysBAM{
template class INCOMPRESSIBLE_EXAMPLE<VECTOR<float,1> >;
template class INCOMPRESSIBLE_EXAMPLE<VECTOR<float,2> >;
template class INCOMPRESSIBLE_EXAMPLE<VECTOR<float,3> >;
template class INCOMPRESSIBLE_EXAMPLE<VECTOR<double,1> >;
template class INCOMPRESSIBLE_EXAMPLE<VECTOR<double,2> >;
template class INCOMPRESSIBLE_EXAMPLE<VECTOR<double,3> >;
}
