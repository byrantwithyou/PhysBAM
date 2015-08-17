//#####################################################################
// Copyright 2009-2010, Michael Lentine, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Parallel_Computation/DOMAIN_ITERATOR_THREADED.h>
#include <Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <Hybrid_Methods/Iterators/PARTICLE_GRID_FACE_WEIGHTS_SPLINE.h>
#include <Hybrid_Methods/Iterators/PARTICLE_GRID_WEIGHTS.h>
#include "SMOKE_EXAMPLE.h"
#include "SMOKE_PARTICLES.h"
#include <pthread.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> SMOKE_EXAMPLE<TV>::
SMOKE_EXAMPLE(const STREAM_TYPE stream_type_input,int number_of_threads)
    :stream_type(stream_type_input),
    debug_particles(*new DEBUG_PARTICLES<TV>),ghost(3),
    initial_time(0),first_frame(0),last_frame(100),frame_rate(24),
    restart(0),write_debug_data(true),output_directory("output"),N_boundary(false),
    debug_divergence(false),alpha(0.1),
    cfl(.9),mac_grid(TV_INT(),RANGE<TV>::Unit_Box(),true),mpi_grid(0),
    thread_queue(number_of_threads>1?new THREAD_QUEUE(number_of_threads):0),projection(mac_grid,false,false,thread_queue),boundary(0),
    use_eapic(false),eapic_order(1),weights(0),particles(*new SMOKE_PARTICLES<TV>)
{
    for(int i=0;i<TV::dimension;i++){domain_boundary(i)(0)=true;domain_boundary(i)(1)=true;}
    pthread_mutex_init(&lock,0);    
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> SMOKE_EXAMPLE<TV>::
~SMOKE_EXAMPLE()
{
    if(mpi_grid || thread_queue) delete boundary;
    delete thread_queue;    
    delete &debug_particles;
    delete weights;
    delete &particles;
}
//#####################################################################
// Function CFL
//#####################################################################
template<class TV> typename TV::SCALAR SMOKE_EXAMPLE<TV>::
CFL(ARRAY<T,FACE_INDEX<TV::dimension> >& face_velocities)
{
    T dt=FLT_MAX;
    DOMAIN_ITERATOR_THREADED_ALPHA<SMOKE_EXAMPLE<TV>,TV>(mac_grid.Domain_Indices(),thread_queue).template Run<ARRAY<T,FACE_INDEX<TV::dimension> >&,T&>(*this,&SMOKE_EXAMPLE::CFL_Threaded,face_velocities,dt);
    return dt;
}
//#####################################################################
// Function CFL_Threaded
//#####################################################################
template<class TV> void SMOKE_EXAMPLE<TV>::
CFL_Threaded(RANGE<TV_INT>& domain,ARRAY<T,FACE_INDEX<TV::dimension> >& face_velocities,T& dt)
{
    T dt_convection=0;
    for(CELL_ITERATOR<TV> iterator(mac_grid,domain);iterator.Valid();iterator.Next()){
        TV_INT cell=iterator.Cell_Index();T local_V_norm=0;
        for(int axis=0;axis<GRID<TV>::dimension;axis++)
            local_V_norm+=mac_grid.one_over_dX[axis]*maxabs(face_velocities(axis,mac_grid.First_Face_Index_In_Cell(axis,cell)),face_velocities(axis,mac_grid.Second_Face_Index_In_Cell(axis,cell)));
        dt_convection=max(dt_convection,local_V_norm);}
    pthread_mutex_lock(&lock);
    dt=min(dt,(T)1.0/dt_convection);
    pthread_mutex_unlock(&lock);
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
    mac_grid.Initialize(counts,domain,true);
}
//#####################################################################
// Function Initialize_Fields
//#####################################################################
template<class TV> void SMOKE_EXAMPLE<TV>::
Initialize_Fields() 
{
    for(FACE_ITERATOR<TV> iterator(mac_grid);iterator.Valid();iterator.Next()) face_velocities(iterator.Full_Index())=0;
    for(CELL_ITERATOR<TV> iterator(mac_grid);iterator.Valid();iterator.Next()) density(iterator.Cell_Index())=0;
}
//#####################################################################
// Function Get_Scalar_Field_Sources
//#####################################################################
template<class TV> void SMOKE_EXAMPLE<TV>::
Get_Scalar_Field_Sources(const T time)
{
    for(CELL_ITERATOR<TV> iterator(mac_grid);iterator.Valid();iterator.Next())
        if(source.Lazy_Inside(iterator.Location())) density(iterator.Cell_Index())=1;
}
//#####################################################################
// Function Set_Weights
//#####################################################################
template<class TV> void SMOKE_EXAMPLE<TV>::
Set_Weights(PARTICLE_GRID_WEIGHTS<TV>* weights_input)
{
    weights=weights_input;
    for(int i=0;i<TV::m;++i){
        if(weights->Order()==1){
            face_weights(i)=new PARTICLE_GRID_FACE_WEIGHTS_SPLINE<TV,1>(mac_grid,/*threads*/1,i);
            face_weights0(i)=new PARTICLE_GRID_FACE_WEIGHTS_SPLINE<TV,1>(mac_grid,/*threads*/1,i);}
        else if(weights->Order()==2){
            face_weights(i)=new PARTICLE_GRID_FACE_WEIGHTS_SPLINE<TV,2>(mac_grid,/*threads*/1,i);
            face_weights0(i)=new PARTICLE_GRID_FACE_WEIGHTS_SPLINE<TV,2>(mac_grid,/*threads*/1,i);}
        else if(weights->Order()==3){
            face_weights(i)=new PARTICLE_GRID_FACE_WEIGHTS_SPLINE<TV,3>(mac_grid,/*threads*/1,i);
            face_weights0(i)=new PARTICLE_GRID_FACE_WEIGHTS_SPLINE<TV,3>(mac_grid,/*threads*/1,i);}
        else PHYSBAM_FATAL_ERROR("Unrecognized interpolation order");}
}
//#####################################################################
// Function Set_Boundary_Conditions
//#####################################################################
template<class TV> void SMOKE_EXAMPLE<TV>::
Set_Boundary_Conditions(const T time)
{
    projection.elliptic_solver->psi_D.Fill(false);projection.elliptic_solver->psi_N.Fill(false);
    for(int axis=0;axis<TV::dimension;axis++) for(int axis_side=0;axis_side<2;axis_side++){int side=2*axis+axis_side;
        if(domain_boundary(axis)(axis_side)){
            TV_INT interior_cell_offset=axis_side==0?TV_INT():-TV_INT::Axis_Vector(axis);    
            for(FACE_ITERATOR<TV> iterator(mac_grid,1,GRID<TV>::BOUNDARY_REGION,side);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Face_Index()+interior_cell_offset;
                TV_INT boundary_face=axis_side==0?iterator.Face_Index()+TV_INT::Axis_Vector(axis):iterator.Face_Index()-TV_INT::Axis_Vector(axis);
                if(N_boundary){
                    FACE_INDEX<TV::dimension> face(axis,boundary_face);
                    projection.elliptic_solver->psi_N(face)=true;
                    face_velocities(face)=0;}
                else projection.elliptic_solver->psi_D(cell)=true;projection.p(cell)=0;}}}
    for(FACE_ITERATOR<TV> iterator(mac_grid);iterator.Valid();iterator.Next()){
        if(source.Lazy_Inside(iterator.Location())){
            projection.elliptic_solver->psi_N(iterator.Full_Index())=true;
            if(iterator.Axis()==1)face_velocities(iterator.Full_Index())=1;
            else face_velocities(iterator.Full_Index())=0;}}
}
//#####################################################################
// Function Write_Output_Files
//#####################################################################
template<class TV> void SMOKE_EXAMPLE<TV>::
Write_Output_Files(const int frame)
{
    std::string f=LOG::sprintf("%d",frame);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/mac_velocities",face_velocities);
    if(mpi_grid) FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/common/global_grid",mpi_grid->global_grid);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/grid",mac_grid);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/common/grid",mac_grid);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/density",density);
    if(write_debug_data){
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/pressure",projection.p);
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/psi_N",projection.elliptic_solver->psi_N);
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/psi_D",projection.elliptic_solver->psi_D);}
    debug_particles.Write_Debug_Particles(stream_type,output_directory,frame);
}
template<class TV> void SMOKE_EXAMPLE<TV>::
Read_Output_Files(const int frame)
{
    std::string f=LOG::sprintf("%d",frame);
    FILE_UTILITIES::Read_From_File(stream_type,output_directory+"/"+f+"/density",density);
    std::string filename;
    filename=output_directory+"/"+f+"/mac_velocities";
    if(FILE_UTILITIES::File_Exists(filename)){LOG::cout<<"Reading mac_velocities "<<filename<<std::endl;FILE_UTILITIES::Read_From_File(stream_type,filename,face_velocities);}
    filename=output_directory+"/"+f+"/pressure";
    if(FILE_UTILITIES::File_Exists(filename)){LOG::cout<<"Reading pressure "<<filename<<std::endl;FILE_UTILITIES::Read_From_File(stream_type,filename,projection.p);}
}
//#####################################################################
namespace PhysBAM{
template class SMOKE_EXAMPLE<VECTOR<float,2> >;
template class SMOKE_EXAMPLE<VECTOR<float,3> >;
template class SMOKE_EXAMPLE<VECTOR<double,2> >;
template class SMOKE_EXAMPLE<VECTOR<double,3> >;
}
