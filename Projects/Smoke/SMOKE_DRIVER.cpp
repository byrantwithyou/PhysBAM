//#####################################################################
// Copyright 2009-2010, Michael Lentine, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Log/DEBUG_SUBSTEPS.h>
#include <Tools/Log/LOG.h>
#include <Tools/Log/SCOPE.h>
#include <Tools/Parallel_Computation/BOUNDARY_MPI.h>
#include <Tools/Parallel_Computation/BOUNDARY_THREADED.h>
#include "SMOKE_DRIVER.h"
#include "SMOKE_EXAMPLE.h"
#include "SMOKE_PARTICLES.h"

using namespace PhysBAM;
namespace{
template<class TV> void Write_Substep_Helper(void* writer,const std::string& title,int substep,int level)
{
    ((SMOKE_DRIVER<TV>*)writer)->Write_Substep(title,substep,level);
}
};
//#####################################################################
// Constructor
//#####################################################################
template<class TV> SMOKE_DRIVER<TV>::
SMOKE_DRIVER(SMOKE_EXAMPLE<TV>& example)
    :ghost(5),example(example)
{
    DEBUG_SUBSTEPS::Set_Substep_Writer((void*)this,&Write_Substep_Helper<TV>);
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> SMOKE_DRIVER<TV>::
~SMOKE_DRIVER()
{
    DEBUG_SUBSTEPS::Clear_Substep_Writer((void*)this);
}
//#####################################################################
// Function Execute_Main_Program
//#####################################################################
template<class TV> void SMOKE_DRIVER<TV>::
Execute_Main_Program()
{
    Initialize();
    Simulate_To_Frame(example.last_frame);
}
//#####################################################################
// Function Initialize
//#####################################################################
template<class TV> void SMOKE_DRIVER<TV>::
Initialize()
{
    DEBUG_SUBSTEPS::Set_Write_Substeps_Level(example.write_substeps_level);

    // setup time
    if(example.restart) current_frame=example.restart;else current_frame=example.first_frame;    
    time=example.Time_At_Frame(current_frame);

    // mpi
    if(example.mpi_grid) example.mpi_grid->Initialize(example.domain_boundary);
    example.projection.elliptic_solver->mpi_grid=example.mpi_grid;
    if(example.mpi_grid) example.boundary=new BOUNDARY_MPI<TV>(example.mpi_grid,example.boundary_scalar);
    else if(example.thread_queue) example.boundary=new BOUNDARY_THREADED<TV>(*example.thread_queue,example.boundary_scalar);    
    else example.boundary=&example.boundary_scalar;

    //threading
    example.projection.elliptic_solver->thread_queue=example.thread_queue;

    // setup grids and velocities
    example.projection.Initialize_Grid(example.mac_grid);
    example.face_velocities.Resize(example.mac_grid);
    example.density.Resize(example.mac_grid.Domain_Indices(ghost));
    example.Initialize_Fields();

    // setup laplace
    example.projection.elliptic_solver->Set_Relative_Tolerance(1e-9);
    example.projection.elliptic_solver->pcg.Set_Maximum_Iterations(1000);
    example.projection.elliptic_solver->pcg.evolution_solver_type=krylov_solver_cg;
    example.projection.elliptic_solver->pcg.cg_restart_iterations=40;

    if(example.restart) example.Read_Output_Files(example.restart);
    
    // setup domain boundaries
    VECTOR<VECTOR<bool,2>,TV::dimension> constant_extrapolation;constant_extrapolation.Fill(VECTOR<bool,2>::Constant_Vector(true));
    example.boundary->Set_Constant_Extrapolation(constant_extrapolation);
    example.Set_Boundary_Conditions(time); // get so CFL is correct

    // EAPIC
    if(example.use_eapic){
        // Resize mass array
        example.mass.Resize(example.mac_grid.Domain_Indices(example.ghost));
        example.face_mass.Resize(example.mac_grid.Domain_Indices(example.ghost));
        // Compute born weights
        for(int i=0;i<TV::m;++i) example.face_weights0(i)->Update(example.particles.X);}

    if(!example.restart) Write_Output_Files(example.first_frame);
    output_number=example.first_frame;
}
//#####################################################################
// Function Update_Particle_Weights
//#####################################################################
template<class TV> void SMOKE_DRIVER<TV>::
Update_Particle_Weights()
{
    PHYSBAM_ASSERT(example.use_eapic);
    example.weights->Update(example.particles.X);
    for(int i=0;i<TV::m;++i) example.face_weights(i)->Update(example.particles.X);
}
//#####################################################################
// Function Particle_To_Grid
//#####################################################################
template<class TV> void SMOKE_DRIVER<TV>::
Particle_To_Grid()
{
    PHYSBAM_ASSERT(example.use_eapic);
    SMOKE_PARTICLES<TV>& particles=example.particles;
    // Zero out grid quantities
    example.mass.Fill((T)0);
    example.face_mass.Fill((T)0);
    example.face_velocities.Fill((T)0);
    // Rasterize mass and momentum to faces
    typename PARTICLE_GRID_ITERATOR<TV>::SCRATCH scratch;
    for(int p=0;p<particles.number;p++){
        for(int d=0;d<TV::m;++d)
            for(PARTICLE_GRID_ITERATOR<TV> it(example.face_weights(d),p,true,scratch);it.Valid();it.Next()){
                T w=it.Weight();
                FACE_INDEX<TV::m> face_index(d,it.Index());
                example.face_mass(face_index)+=w*particles.mass(p);
                TV dx=example.mac_grid.Face(face_index)-particles.X(p);
                T v=particles.V(p)(d)+TV::Dot_Product(particles.C(p).Column(d),dx);
                example.face_velocities(face_index)+=w*particles.mass(p)*v;}}
    // Divide out masses
    for(FACE_ITERATOR<TV> iterator(example.mac_grid,ghost);iterator.Valid();iterator.Next()){
        FACE_INDEX<TV::m> face_index=iterator.Full_Index();
        if(example.face_mass(face_index))
            example.face_velocities(face_index)/=example.face_mass(face_index);}
    // TODO: Cell center stuff (for density/temperature advection)
}
//#####################################################################
// Function Grid_To_Particle
//#####################################################################
template<class TV> void SMOKE_DRIVER<TV>::
Grid_To_Particle()
{
    SMOKE_PARTICLES<TV>& particles=example.particles;
    typename PARTICLE_GRID_ITERATOR<TV>::SCRATCH scratch;

#pragma omp parallel for
    for(int p=0;p<particles.number;p++){
        // Zero out particle veloicity and C, resample particle positions
        particles.X(p)=particles.X0(p);
        particles.V(p)=TV();
        particles.C(p)=MATRIX<T,TV::m>();
        // Compute new V and C
        for(int d=0;d<TV::m;++d)
            for(PARTICLE_GRID_ITERATOR<TV> it(example.face_weights(d),p,true,scratch);it.Valid();it.Next()){
                T w=it.Weight();
                FACE_INDEX<TV::m> face_index(d,it.Index());
                particles.V(p)(d)+=w*example.face_velocities(face_index);
                particles.C(p).Add_Column(d,example.face_velocities(face_index)*it.Gradient());}}
}
//#####################################################################
// Function Move_Particles
//#####################################################################
template<class TV> void SMOKE_DRIVER<TV>::
Move_Particles(const T dt)
{
    PHYSBAM_ASSERT(example.use_eapic);
    SMOKE_PARTICLES<TV>& particles=example.particles;
#pragma omp parallel for
    for(int p=0;p<particles.number;p++)
        particles.X(p)+=particles.V(p)*dt;
}
//#####################################################################
// Function Add_Buoyancy_Force
//#####################################################################
template<class TV> void SMOKE_DRIVER<TV>::
Add_Buoyancy_Force(const T dt,const T time)
{
    for(FACE_ITERATOR<TV> iterator(example.mac_grid);iterator.Valid();iterator.Next()){
        if(iterator.axis==1){
            TV_INT c1,c2;
            example.mac_grid.Cells_Touching_Face(iterator.axis,iterator.Face_Index(),c1,c2);
            T rho=(example.density(c1)+example.density(c2))*(T).5;
            // T tem=(example.temperature(c1)+example.temperature(c2))*(T).5; // no temperature for now 
            example.face_velocities(iterator.Full_Index())+=-dt*example.alpha*rho;}}
}
//#####################################################################
// Function Scalar_Advance
//#####################################################################
template<class TV> void SMOKE_DRIVER<TV>::
Scalar_Advance(const T dt,const T time)
{
    example.Get_Scalar_Field_Sources(time);
    ARRAY<T,TV_INT> density_ghost(example.mac_grid.Domain_Indices(ghost));
    example.boundary->Set_Fixed_Boundary(true,0);
    example.boundary->Fill_Ghost_Cells(example.mac_grid,example.density,density_ghost,dt,time,ghost);
    example.advection_scalar.Update_Advection_Equation_Cell(example.mac_grid,example.density,density_ghost,example.face_velocities,*example.boundary,dt,time);    
    example.boundary->Set_Fixed_Boundary(false);
}
//#####################################################################
// Function Convect
//#####################################################################
template<class TV> void SMOKE_DRIVER<TV>::
Convect(const T dt,const T time)
{
    if(example.use_eapic){
        Grid_To_Particle();
        Move_Particles(dt);
        Update_Particle_Weights();
        Particle_To_Grid();}
    else{
        example.boundary->Set_Fixed_Boundary(true,0);
        ARRAY<T,FACE_INDEX<TV::dimension> > face_velocities_ghost(example.mac_grid,ghost,false);
        example.boundary->Fill_Ghost_Faces(example.mac_grid,example.face_velocities,face_velocities_ghost,time,ghost);
        example.advection_scalar.Update_Advection_Equation_Face(
            example.mac_grid,example.face_velocities,face_velocities_ghost,face_velocities_ghost,*example.boundary,dt,time);
        example.boundary->Set_Fixed_Boundary(false);}
}
//#####################################################################
// Function Project
//#####################################################################
template<class TV> void SMOKE_DRIVER<TV>::
Project(const T dt,const T time)
{
    example.Set_Boundary_Conditions(time+dt);
    example.projection.p*=dt; // rescale pressure for guess
    example.boundary->Apply_Boundary_Condition_Face(example.mac_grid,example.face_velocities,time+dt);
    example.projection.Make_Divergence_Free(example.face_velocities,dt,time);
    example.projection.p*=(1/dt); // unscale pressure
}
//#####################################################################
// Function Advance_To_Target_Time
//#####################################################################
template<class TV> void SMOKE_DRIVER<TV>::
Advance_To_Target_Time(const T target_time)
{
    bool done=false;for(int substep=1;!done;substep++){
        // choose time step
        LOG::Time("Substep");
        T dt=example.cfl*example.CFL(example.face_velocities);        
        if(example.mpi_grid) example.mpi_grid->Synchronize_Dt(dt);
        if(time+dt>=target_time){dt=target_time-time;done=true;}
        else if(time+2*dt>=target_time){dt=.5*(target_time-time);}
        Scalar_Advance(dt,time);
        Convect(dt,time);
        Add_Buoyancy_Force(dt,time);
        Print_Max_Divergence("Before projection");
        Project(dt,time);
        Print_Max_Divergence("After projection");
        time+=dt;}
}
//#####################################################################
// Function Simulate_To_Frame
//#####################################################################
template<class TV> void SMOKE_DRIVER<TV>::
Simulate_To_Frame(const int frame)
{
    while(current_frame<frame){
        LOG::SCOPE scope("FRAME","Frame %d",current_frame+1);
        Advance_To_Target_Time(example.Time_At_Frame(current_frame+1));
        Write_Output_Files(++output_number);
        current_frame++;}
}
//#####################################################################
// Function Write_Substep
//#####################################################################
template<class TV> void SMOKE_DRIVER<TV>::
Write_Substep(const std::string& title,const int substep,const int level)
{
    if(level<=example.write_substeps_level){
        example.frame_title=title;
        LOG::cout<<"Writing substep ["<<example.frame_title<<"]: output_number="<<output_number+1<<", time="<<time<<", frame="<<current_frame<<", substep="<<substep<<std::endl;
        Write_Output_Files(++output_number);example.frame_title="";}
}
//#####################################################################
// Function Write_Output_Files
//#####################################################################
template<class TV> void SMOKE_DRIVER<TV>::
Write_Output_Files(const int frame)
{
    FILE_UTILITIES::Create_Directory(example.output_directory);
    FILE_UTILITIES::Create_Directory(example.output_directory+LOG::sprintf("/%d",frame));
    FILE_UTILITIES::Create_Directory(example.output_directory+"/common");
    FILE_UTILITIES::Write_To_Text_File(example.output_directory+LOG::sprintf("/%d/frame_title",frame),example.frame_title);
    if(frame==example.first_frame) 
        FILE_UTILITIES::Write_To_Text_File(example.output_directory+"/common/first_frame",frame,"\n");
    example.Write_Output_Files(frame);
    FILE_UTILITIES::Write_To_Text_File(example.output_directory+"/common/last_frame",frame,"\n");
}
//#####################################################################
// Function Print_Max_Divergence
//#####################################################################
template<class TV> void SMOKE_DRIVER<TV>::
Print_Max_Divergence(const char* str)
{
    if(!example.debug_divergence) return;
    T max_div=(T)0;
    for(CELL_ITERATOR<TV> iterator(example.mac_grid);iterator.Valid();iterator.Next()){
        T d=0;
        for(int axis=0;axis<TV::m;axis++)
            d+=example.face_velocities(FACE_INDEX<TV::m>(axis,iterator.Second_Face_Index(axis)))
                -example.face_velocities(FACE_INDEX<TV::m>(axis,iterator.First_Face_Index(axis)));
        T ad=abs(d); if(ad>max_div) max_div=ad;}
    LOG::cout<<str<<" max(div(v)) "<<max_div<<std::endl;
}
//#####################################################################
namespace PhysBAM{
template class SMOKE_DRIVER<VECTOR<float,2> >;
template class SMOKE_DRIVER<VECTOR<float,3> >;
template class SMOKE_DRIVER<VECTOR<double,2> >;
template class SMOKE_DRIVER<VECTOR<double,3> >;
}
