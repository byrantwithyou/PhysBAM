//#####################################################################
// Copyright 2009-2010, Michael Lentine, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Log/LOG.h>
#include <Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <Hybrid_Methods/Iterators/PARTICLE_GRID_WEIGHTS.h>
#include <Hybrid_Methods/Iterators/PARTICLE_GRID_WEIGHTS_SPLINE.h>
#include "SMOKE_EXAMPLE.h"
#include "SMOKE_PARTICLES.h"
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> SMOKE_EXAMPLE<TV>::
SMOKE_EXAMPLE(const STREAM_TYPE stream_type_input,int number_of_threads)
    :stream_type(stream_type_input),
    debug_particles(*new DEBUG_PARTICLES<TV>),ghost(5),
    initial_time(0),first_frame(0),last_frame(100),frame_rate(24),
    restart(0),write_debug_data(true),output_directory("output"),N_boundary(false),
    debug_divergence(false),alpha(0.1),beta(0.00366),
    cfl(.9),grid(TV_INT(),RANGE<TV>::Unit_Box(),true),mpi_grid(0),
    projection(grid,false,false),boundary(0),
    use_eapic(false),eapic_order(1),particles(*new SMOKE_PARTICLES<TV>)
{
    np=1; //number of points per cell
    for(int i=0;i<TV::m;i++){domain_boundary(i)(0)=true;domain_boundary(i)(1)=true;}
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> SMOKE_EXAMPLE<TV>::
~SMOKE_EXAMPLE()
{
    if(mpi_grid) delete boundary;
    delete &debug_particles;
    for(int i=0;i<TV::m;i++) delete weights(i);
    delete &particles;
}
//#####################################################################
// Function CFL
//#####################################################################
template<class TV> typename TV::SCALAR SMOKE_EXAMPLE<TV>::
CFL(ARRAY<T,FACE_INDEX<TV::m> >& face_velocities)
{
    T dt=FLT_MAX;
    T dt_convection=0;
    for(CELL_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()){
        TV_INT cell=iterator.Cell_Index();T local_V_norm=0;
        for(int axis=0;axis<GRID<TV>::dimension;axis++)
            local_V_norm+=grid.one_over_dX[axis]*maxabs(face_velocities(axis,grid.First_Face_Index_In_Cell(axis,cell)),face_velocities(axis,grid.Second_Face_Index_In_Cell(axis,cell)));
        dt_convection=max(dt_convection,local_V_norm);}
    dt=min(dt,(T)1.0/dt_convection);
    return dt;
}
//#####################################################################
// Function Time_At_Frame
//#####################################################################
template<class TV> typename TV::SCALAR SMOKE_EXAMPLE<TV>::
Time_At_Frame(const int frame) const 
{
    return initial_time+(frame-first_frame)/frame_rate;
}
//#####################################################################
// Function Initialize_Grid
//#####################################################################
template<class TV> void SMOKE_EXAMPLE<TV>::
Initialize_Grid(TV_INT counts,RANGE<TV> domain) 
{
    grid.Initialize(counts,domain,true);
}
//#####################################################################
// Function Initialize_Fields
//#####################################################################
template<class TV> void SMOKE_EXAMPLE<TV>::
Initialize_Fields() 
{
    for(FACE_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()) face_velocities(iterator.Full_Index())=0;
    for(CELL_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()) density(iterator.Cell_Index())=0;
    for(CELL_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()) temperature(iterator.Cell_Index())=273; // add temperature
}
//#####################################################################
// Function Get_Scalar_Field_Sources
//#####################################################################
template<class TV> void SMOKE_EXAMPLE<TV>::
Get_Scalar_Field_Sources(const T time)
{
    for(CELL_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next())
    {
        if(source1.Lazy_Inside(iterator.Location())) {density(iterator.Cell_Index())=1;}
        if(source1.Lazy_Inside(iterator.Location())) {temperature(iterator.Cell_Index())=275;}    
        //else if(source2.Lazy_Inside(iterator.Location())) {density(iterator.Cell_Index())=1;}
    }
}
//#####################################################################
// Function Set_Weights
//#####################################################################
template<class TV> void SMOKE_EXAMPLE<TV>::
Set_Weights(int order)
{
    for(int i=0;i<TV::m;++i){
        GRID<TV> face_grid=grid.Get_Face_Grid(i).Get_MAC_Grid_At_Regular_Positions();
        if(order==1){
            weights(i)=new PARTICLE_GRID_WEIGHTS_SPLINE<TV,1>(grid,/*threads*/1);
            weights0(i)=new PARTICLE_GRID_WEIGHTS_SPLINE<TV,1>(grid,/*threads*/1);}
        else if(order==2){
            weights(i)=new PARTICLE_GRID_WEIGHTS_SPLINE<TV,2>(grid,/*threads*/1);
            weights0(i)=new PARTICLE_GRID_WEIGHTS_SPLINE<TV,2>(grid,/*threads*/1);}
        else if(order==3){
            weights(i)=new PARTICLE_GRID_WEIGHTS_SPLINE<TV,3>(grid,/*threads*/1);
            weights0(i)=new PARTICLE_GRID_WEIGHTS_SPLINE<TV,3>(grid,/*threads*/1);}
        else PHYSBAM_FATAL_ERROR("Unrecognized interpolation order");
        weights(i)->use_gradient_transfer=true;
        weights0(i)->use_gradient_transfer=true;}
}
//#####################################################################
// Function Set_Boundary_Conditions
//#####################################################################
template<class TV> void SMOKE_EXAMPLE<TV>::
Set_Boundary_Conditions(const T time, const T source_velocities)
{
    projection.elliptic_solver->psi_D.Fill(false);projection.elliptic_solver->psi_N.Fill(false);
    for(int axis=0;axis<TV::m;axis++) for(int axis_side=0;axis_side<2;axis_side++){int side=2*axis+axis_side;
        if(domain_boundary(axis)(axis_side)){
            TV_INT interior_cell_offset=axis_side==0?TV_INT():-TV_INT::Axis_Vector(axis);    
            for(FACE_ITERATOR<TV> iterator(grid,1,GRID<TV>::BOUNDARY_REGION,side);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Face_Index()+interior_cell_offset;
                TV_INT boundary_face=axis_side==0?iterator.Face_Index()+TV_INT::Axis_Vector(axis):iterator.Face_Index()-TV_INT::Axis_Vector(axis);
                if(N_boundary){
                    FACE_INDEX<TV::m> face(axis,boundary_face);
                    projection.elliptic_solver->psi_N(face)=true;
                    face_velocities(face)=0;}
                else projection.elliptic_solver->psi_D(cell)=true;
                projection.p(cell)=0;}}}
    for(FACE_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()){
        if(source1.Lazy_Inside(iterator.Location())){
            projection.elliptic_solver->psi_N(iterator.Full_Index())=true;
            if(iterator.Axis()==1)face_velocities(iterator.Full_Index())=source_velocities;
            else face_velocities(iterator.Full_Index())=0;}
        // else if(source2.Lazy_Inside(iterator.Location())){   // second source
        //     projection.elliptic_solver->psi_N(iterator.Full_Index())=true;
        //     if(iterator.Axis()==1)face_velocities(iterator.Full_Index())=-source_velocities;
        //     else face_velocities(iterator.Full_Index())=0;}
    }
}
//#####################################################################
// Function Write_Output_Files
//#####################################################################
template<class TV> void SMOKE_EXAMPLE<TV>::
Write_Output_Files(const int frame)
{
    std::string f=LOG::sprintf("%d",frame);
    Write_To_File(stream_type,output_directory+"/"+f+"/mac_velocities",face_velocities);
    if(mpi_grid) Write_To_File(stream_type,output_directory+"/common/global_grid",mpi_grid->global_grid);
    Write_To_File(stream_type,output_directory+"/"+f+"/grid",grid);
    Write_To_File(stream_type,output_directory+"/common/grid",grid);
    Write_To_File(stream_type,output_directory+"/"+f+"/density",density);
    Write_To_File(stream_type,output_directory+"/"+f+"/temperature",temperature);// add temperature
    if(write_debug_data){
        Write_To_File(stream_type,output_directory+"/"+f+"/pressure",projection.p);
        Write_To_File(stream_type,output_directory+"/"+f+"/psi_N",projection.elliptic_solver->psi_N);
        Write_To_File(stream_type,output_directory+"/"+f+"/psi_D",projection.elliptic_solver->psi_D);}
    for(int p=0;p<particles.number;p++){
        Add_Debug_Particle(particles.X(p),VECTOR<T,3>(0,1,0));
        Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_V,particles.V(p));}
    debug_particles.Write_Debug_Particles(stream_type,output_directory,frame);
}
template<class TV> void SMOKE_EXAMPLE<TV>::
Read_Output_Files(const int frame)
{
    std::string f=LOG::sprintf("%d",frame);
    Read_From_File(stream_type,output_directory+"/"+f+"/density",density);
    std::string filename;
    filename=output_directory+"/"+f+"/mac_velocities";
    if(File_Exists(filename)){LOG::cout<<"Reading mac_velocities "<<filename<<std::endl;Read_From_File(stream_type,filename,face_velocities);}
    filename=output_directory+"/"+f+"/pressure";
    if(File_Exists(filename)){LOG::cout<<"Reading pressure "<<filename<<std::endl;Read_From_File(stream_type,filename,projection.p);}
}
//#####################################################################
namespace PhysBAM{
template class SMOKE_EXAMPLE<VECTOR<float,2> >;
template class SMOKE_EXAMPLE<VECTOR<float,3> >;
template class SMOKE_EXAMPLE<VECTOR<double,2> >;
template class SMOKE_EXAMPLE<VECTOR<double,3> >;
}
