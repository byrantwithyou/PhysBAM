//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Arrays/INDIRECT_ARRAY.h>
#include <Tools/Grids_Uniform/FACE_ITERATOR.h>
#include <Tools/Grids_Uniform/NODE_ITERATOR.h>
#include <Tools/Krylov_Solvers/GMRES.h>
#include <Tools/Krylov_Solvers/KRYLOV_SOLVER.h>
#include <Tools/Krylov_Solvers/MINRES.h>
#include <Tools/Log/DEBUG_SUBSTEPS.h>
#include <Tools/Log/LOG.h>
#include <Tools/Nonlinear_Equations/NEWTONS_METHOD.h>
#include <Tools/Read_Write/OCTAVE_OUTPUT.h>
#include <Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_DRIVER.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_EXAMPLE.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_PARTICLES.h>
#include <Hybrid_Methods/Iterators/GATHER_SCATTER.h>
#include <Hybrid_Methods/Iterators/PARTICLE_GRID_WEIGHTS.h>
#include <Hybrid_Methods/System/MPM_KRYLOV_VECTOR.h>
#include <Hybrid_Methods/System/MPM_OBJECTIVE.h>
#include <boost/function.hpp>
using namespace PhysBAM;
namespace{
    template<class TV> void Write_Substep_Helper(void* writer,const std::string& title,int substep,int level)
    {
        ((MPM_DRIVER<TV>*)writer)->Write_Substep(title,substep,level);
    }
};
//#####################################################################
// Constructor
//#####################################################################
template<class TV> MPM_DRIVER<TV>::
MPM_DRIVER(MPM_EXAMPLE<TV>& example)
    :example(example),objective(*new MPM_OBJECTIVE<TV>(example)),
    dv(*new MPM_KRYLOV_VECTOR<TV>(example.valid_grid_indices)),
    rhs(*new MPM_KRYLOV_VECTOR<TV>(example.valid_grid_indices))
{
    DEBUG_SUBSTEPS::Set_Substep_Writer((void*)this,&Write_Substep_Helper<TV>);
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> MPM_DRIVER<TV>::
~MPM_DRIVER()
{
    DEBUG_SUBSTEPS::Clear_Substep_Writer((void*)this);
    delete &objective;
}
//#####################################################################
// Execute_Main_Program
//#####################################################################
template<class TV> void MPM_DRIVER<TV>::
Execute_Main_Program()
{
    Initialize();
    Simulate_To_Frame(example.last_frame);
}
//#####################################################################
// Initialize
//#####################################################################
template<class TV> void MPM_DRIVER<TV>::
Initialize()
{
    LOG::cout<<std::setprecision(16)<<std::endl;
    DEBUG_SUBSTEPS::Set_Write_Substeps_Level(example.substeps_delay_frame<0?example.write_substeps_level:-1);

    // setup time
    output_number=current_frame=example.restart;

    example.Initialize();
    if(example.restart)
        example.Read_Output_Files(example.restart);

    example.mass.Resize(example.grid.Domain_Indices(example.ghost));
    example.velocity.Resize(example.grid.Domain_Indices(example.ghost));
    example.rhs.u.Resize(example.grid.Domain_Indices(example.ghost));

    if(!example.restart) Write_Output_Files(0);
    PHYSBAM_DEBUG_WRITE_SUBSTEP("after init",0,1);
}
//#####################################################################
// Function Advance_One_Time_Step
//#####################################################################
template<class TV> void MPM_DRIVER<TV>::
Advance_One_Time_Step()
{
    example.Begin_Time_Step(example.time);

    Update_Particle_Weights();

    Particle_To_Grid();

    Apply_Forces();

    Grid_To_Particle();

    example.End_Time_Step(example.time);
}
// Simulate_To_Frame
//#####################################################################
template<class TV> void MPM_DRIVER<TV>::
Simulate_To_Frame(const int frame)
{
    for(;current_frame<frame;current_frame++){
        LOG::SCOPE scope("FRAME","frame %d",current_frame+1);
        example.Begin_Frame(current_frame);
        if(example.substeps_delay_frame==current_frame)
            DEBUG_SUBSTEPS::Set_Write_Substeps_Level(example.write_substeps_level);
        T time_at_frame=example.time+example.frame_dt;
        bool done=false;
        for(int substep=0;!done;substep++){
            LOG::SCOPE scope("SUBSTEP","substep %d",substep+1);
            example.dt=Compute_Dt();
            // example.dt=example.max_dt;

            example.dt=clamp(example.dt,example.min_dt,example.max_dt);
            T next_time=example.time+example.dt;
            if(next_time>time_at_frame){
                next_time=time_at_frame;
                done=true;}
            else if(next_time+example.dt>time_at_frame) next_time=(example.time+time_at_frame)/2;
            example.dt=next_time-example.time;

            Advance_One_Time_Step();
            example.time=next_time;}
        example.End_Frame(current_frame);
        Write_Output_Files(++output_number);}
}
//#####################################################################
// Function Write_Substep
//#####################################################################
template<class TV> void MPM_DRIVER<TV>::
Write_Substep(const std::string& title,const int substep,const int level)
{
    if(level<=example.write_substeps_level){
        example.frame_title=title;
        LOG::printf("Writing substep [%s]: output_number=%i, time=%g, frame=%i, substep=%i\n",
            example.frame_title,output_number+1,time,current_frame,substep);
        Write_Output_Files(++output_number);
        example.frame_title="";}
}
//#####################################################################
// Write_Output_Files
//#####################################################################
template<class TV> void MPM_DRIVER<TV>::
Write_Output_Files(const int frame)
{
    FILE_UTILITIES::Create_Directory(example.output_directory);
    FILE_UTILITIES::Create_Directory(example.output_directory+STRING_UTILITIES::string_sprintf("/%d",frame));
    FILE_UTILITIES::Create_Directory(example.output_directory+"/common");
    FILE_UTILITIES::Write_To_Text_File(example.output_directory+STRING_UTILITIES::string_sprintf("/%d/frame_title",frame),example.frame_title);
    if(frame==0)
        FILE_UTILITIES::Write_To_Text_File(example.output_directory+"/common/first_frame",frame,"\n");
    example.Write_Output_Files(frame);
    FILE_UTILITIES::Write_To_Text_File(example.output_directory+"/common/last_frame",frame,"\n");
}
//#####################################################################
// Function Update_Particle_Weights
//#####################################################################
template<class TV> void MPM_DRIVER<TV>::
Update_Particle_Weights()
{
    example.weights->Update(example.particles.X);
}
//#####################################################################
// Function Particle_To_Grid
//#####################################################################
template<class TV> void MPM_DRIVER<TV>::
Particle_To_Grid()
{
    MPM_PARTICLES<TV>& particles=example.particles;

    if(example.weights->use_gradient_transfer)
    {
        example.gather_scatter.Scatter(
            [this,&particles](int p,const PARTICLE_GRID_ITERATOR<TV>& it,int tid)
            {
                example.mass(it.Index())+=it.Weight()*particles.mass(p);
                example.velocity(it.Index())+=particles.mass(p)*(it.Weight()*particles.V(p)+particles.B(p)*it.Gradient());
            },true);
    }
    else if(example.weights->constant_scalar_inertia_tensor)
    {
        T Dp_inverse=example.weights->Constant_Scalar_Inverse_Dp();
        example.gather_scatter.Scatter(
            [this,Dp_inverse,&particles](int p,const PARTICLE_GRID_ITERATOR<TV>& it,int tid)
            {
                example.mass(it.Index())+=it.Weight()*particles.mass(p);
                TV V=particles.V(p);
                if(example.use_affine) V+=particles.B(p)*(Dp_inverse*(example.grid.Center(it.Index())-particles.X(p)));
                example.velocity(it.Index())+=it.Weight()*particles.mass(p)*V;
            },false);
    }
    else PHYSBAM_FATAL_ERROR("General case for rasterization not implemented");

    for(int i=0;i<example.mass.array.m;i++)
        if(example.mass.array(i))
            example.velocity.array(i)/=example.mass.array(i);
        else example.velocity.array(i)=TV();
}
//#####################################################################
// Function Grid_To_Particle
//#####################################################################
template<class TV> void MPM_DRIVER<TV>::
Grid_To_Particle()
{
    MPM_PARTICLES<TV>& particles=example.particles;
    T dt=example.dt;

    for(int k=0;k<example.simulated_particles.m;k++){
        int p=example.simulated_particles(k);
        TV Vn_interpolate,V_pic,V_flip=particles.V(p);
        MATRIX<T,TV::m> B,grad_Vp;

        for(PARTICLE_GRID_ITERATOR<TV> it(example.weights,p,true,0);it.Valid();it.Next()){
            T w=it.Weight();
            TV_INT index=it.Index();
            TV V_grid=example.velocity_new(index);
            V_pic+=w*V_grid;
            V_flip+=w*(V_grid-example.velocity(index));
            Vn_interpolate+=w*example.velocity(index);
            if(example.use_midpoint)
                V_grid=(T).5*(V_grid+example.velocity(index));
            grad_Vp+=MATRIX<T,TV::m>::Outer_Product(V_grid,it.Gradient());}
        particles.F(p)+=dt*grad_Vp*particles.F(p);

        if(example.use_affine)
            for(PARTICLE_GRID_ITERATOR<TV> it(example.weights,p,false,0);it.Valid();it.Next()){
                TV_INT index=it.Index();
                TV V_grid=example.velocity_new(index);
                TV Z=example.grid.Center(index);
                TV xi_new=Z+dt/2*(V_grid+example.velocity(index));
                TV xp_new=particles.X(p)+dt/2*(Vn_interpolate+V_pic);
                B+=it.Weight()/2*(MATRIX<T,TV::m>::Outer_Product(V_grid,Z-particles.X(p)+xi_new-xp_new)
                    +MATRIX<T,TV::m>::Outer_Product(Z-particles.X(p)-xi_new+xp_new,V_grid));}

        particles.V(p)=V_pic;
        Perform_Particle_Collision(p);
        if(example.use_midpoint) particles.X(p)+=(particles.V(p)+Vn_interpolate)*(dt/2);
        else particles.X(p)+=particles.V(p)*dt;
        particles.V(p)=V_flip*example.flip+V_pic*(1-example.flip);
        particles.B(p)=B;}
}
//#####################################################################
// Function Apply_Forces
//#####################################################################
template<class TV> void MPM_DRIVER<TV>::
Apply_Forces()
{
    objective.Reset();
    NEWTONS_METHOD<T> newtons_method;
    newtons_method.tolerance=example.newton_tolerance*example.dt;
    newtons_method.progress_tolerance=1e-5;
    newtons_method.max_iterations=example.newton_iterations;
    newtons_method.krylov_tolerance=example.solver_tolerance;
    newtons_method.max_krylov_iterations=example.solver_iterations;
    newtons_method.use_cg=true;

    newtons_method.require_one_iteration=!objective.Initial_Guess(dv,newtons_method.tolerance);
    if(example.test_diff) objective.Test_Diff(dv);

    bool converged=newtons_method.Newtons_Method(objective,objective.system,dv,av);
    if(!converged) LOG::cout<<"WARNING: Newton's method did not converge"<<std::endl;

    Apply_Friction();
    objective.Restore_F();

    for(int i=0;i<example.valid_grid_indices.m;i++){
        int j=example.valid_grid_indices(i);
        example.velocity_new.array(j)=dv.u.array(j)+objective.v0.u.array(j);}
}
//#####################################################################
// Function Perform_Particle_Collision
//#####################################################################
template<class TV> void MPM_DRIVER<TV>::
Perform_Particle_Collision(int p)
{
    // TODO
}
//#####################################################################
// Function Apply_Friction
//#####################################################################
template<class TV> void MPM_DRIVER<TV>::
Apply_Friction()
{
    // TODO
}
//#####################################################################
// Function Compute_Dt
//#####################################################################
template<class TV> typename TV::SCALAR MPM_DRIVER<TV>::
Compute_Dt() const
{
    T critical_speed=example.cfl*example.grid.One_Over_DX().Min()/example.max_dt;
    T v=Max_Particle_Speed();
    return (v>critical_speed)?(example.cfl*example.grid.One_Over_DX().Min()/v):example.max_dt;
}
//#####################################################################
// Function Max_Particle_Speed
//#####################################################################
template<class TV> typename TV::SCALAR MPM_DRIVER<TV>::
Max_Particle_Speed() const
{
    T v2=0;
    for(int k=0;k<example.simulated_particles.m;k++){
        int p=example.simulated_particles(k);
        v2=max(v2,example.particles.V(p).Magnitude_Squared());}
    return sqrt(v2);
}
//#####################################################################
namespace PhysBAM{
template class MPM_DRIVER<VECTOR<float,1> >;
template class MPM_DRIVER<VECTOR<float,2> >;
template class MPM_DRIVER<VECTOR<float,3> >;
template class MPM_DRIVER<VECTOR<double,1> >;
template class MPM_DRIVER<VECTOR<double,2> >;
template class MPM_DRIVER<VECTOR<double,3> >;
}
