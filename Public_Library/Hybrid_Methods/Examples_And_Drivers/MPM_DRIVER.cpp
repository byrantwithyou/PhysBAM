//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Arrays/INDIRECT_ARRAY.h>
#include <Tools/Grids_Uniform/CELL_ITERATOR.h>
#include <Tools/Grids_Uniform/FACE_ITERATOR.h>
#include <Tools/Grids_Uniform/NODE_ITERATOR.h>
#include <Tools/Krylov_Solvers/GMRES.h>
#include <Tools/Krylov_Solvers/KRYLOV_SOLVER.h>
#include <Tools/Krylov_Solvers/MINRES.h>
#include <Tools/Log/DEBUG_SUBSTEPS.h>
#include <Tools/Log/LOG.h>
#include <Tools/Log/SCOPE.h>
#include <Tools/Nonlinear_Equations/NEWTONS_METHOD.h>
#include <Tools/Read_Write/OCTAVE_OUTPUT.h>
#include <Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>
#include <Deformables/Forces/DEFORMABLES_FORCES.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_DRIVER.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_EXAMPLE.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_PARTICLES.h>
#include <Hybrid_Methods/Iterators/GATHER_SCATTER.h>
#include <Hybrid_Methods/Iterators/PARTICLE_GRID_WEIGHTS.h>
#include <Hybrid_Methods/System/MPM_KRYLOV_VECTOR.h>
#include <Hybrid_Methods/System/MPM_OBJECTIVE.h>
#include <boost/function.hpp>
#ifdef USE_OPENMP
#include <omp.h>
#endif
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
    delete &dv;
    delete &rhs;
    av.Delete_Pointers_And_Clean_Memory();
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
    PHYSBAM_ASSERT(example.grid.Is_MAC_Grid());
    if(example.restart)
        example.Read_Output_Files(example.restart);

    example.mass.Resize(example.grid.Domain_Indices(example.ghost));
    example.velocity.Resize(example.grid.Domain_Indices(example.ghost));
    example.velocity_new.Resize(example.grid.Domain_Indices(example.ghost));
    dv.u.Resize(example.grid.Domain_Indices(example.ghost));
    rhs.u.Resize(example.grid.Domain_Indices(example.ghost));
    objective.system.tmp.u.Resize(example.grid.Domain_Indices(example.ghost));

    RANGE<TV_INT> range(example.grid.Cell_Indices(example.ghost));
    example.location.Resize(range,false,false);
#pragma omp parallel for
    for(int t=0;t<example.threads;t++){
        int a=(range.max_corner.x-range.min_corner.x)*t/example.threads+range.min_corner.x;
        int b=(range.max_corner.x-range.min_corner.x)*(t+1)/example.threads+range.min_corner.x;
        RANGE<TV_INT> local_range(range);
        local_range.max_corner.x=a;
        int i=local_range.Size();
        local_range.min_corner.x=a;
        local_range.max_corner.x=b;
        for(RANGE_ITERATOR<TV::m> it(local_range);it.Valid();it.Next())
            example.location.array(i++)=example.grid.Center(it.index);}

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

    Print_Particle_Stats("particle state",example.dt,example.velocity,0);
    Update_Simulated_Particles();
    Update_Particle_Weights();
    example.gather_scatter.Prepare_Scatter(example.particles);
    Particle_To_Grid();
    Print_Grid_Stats("after particle to grid",example.dt,example.velocity,0);
    Print_Energy_Stats("after particle to grid",example.velocity);
    PHYSBAM_DEBUG_WRITE_SUBSTEP("after particle to grid",0,1);
    Apply_Forces();
    Print_Grid_Stats("after forces",example.dt,example.velocity_new,&example.velocity);
    PHYSBAM_DEBUG_WRITE_SUBSTEP("after forces",0,1);
    Grid_To_Particle();
    PHYSBAM_DEBUG_WRITE_SUBSTEP("after grid to particle",0,1);

    example.End_Time_Step(example.time);
}
//#####################################################################
// Function Simulate_To_Frame
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
            example.dt=clamp(example.dt,example.min_dt,example.max_dt);
            LOG::cout<<"substep dt: "<<example.dt<<std::endl;
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
    FILE_UTILITIES::Create_Directory(example.output_directory+LOG::sprintf("/%d",frame));
    FILE_UTILITIES::Create_Directory(example.output_directory+"/common");
    FILE_UTILITIES::Write_To_Text_File(example.output_directory+LOG::sprintf("/%d/frame_title",frame),example.frame_title);
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

    example.mass.array.Fill(0);
    example.velocity.array.Fill(TV());
    if(example.weights->use_gradient_transfer)
    {
        example.gather_scatter.template Scatter<int>(
            [this,&particles](int p,const PARTICLE_GRID_ITERATOR<TV>& it,int data)
            {
                example.mass(it.Index())+=it.Weight()*particles.mass(p);
                example.velocity(it.Index())+=particles.mass(p)*(it.Weight()*particles.V(p)+particles.B(p)*it.Gradient());
            },true);
    }
    else if(example.weights->constant_scalar_inertia_tensor)
    {
        T Dp_inverse=example.weights->Constant_Scalar_Inverse_Dp();
        example.gather_scatter.template Scatter<int>(
            [this,Dp_inverse,&particles](int p,const PARTICLE_GRID_ITERATOR<TV>& it,int data)
            {
                example.mass(it.Index())+=it.Weight()*particles.mass(p);
                TV V=particles.V(p);
                if(example.use_affine) V+=particles.B(p)*(Dp_inverse*(example.grid.Center(it.Index())-particles.X(p)));
                example.velocity(it.Index())+=it.Weight()*particles.mass(p)*V;
            },false);
    }
    else PHYSBAM_FATAL_ERROR("General case for rasterization not implemented");

    example.valid_grid_indices.Remove_All();
    for(int i=0;i<example.mass.array.m;i++)
        if(example.mass.array(i)){
            example.valid_grid_indices.Append(i);
            example.velocity.array(i)/=example.mass.array(i);}
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

#pragma omp parallel for
    for(int tid=0;tid<example.gather_scatter.threads;++tid){
        int a=example.simulated_particles.m*tid/example.gather_scatter.threads;
        int b=example.simulated_particles.m*(tid+1)/example.gather_scatter.threads;
        typename PARTICLE_GRID_ITERATOR<TV>::SCRATCH scratch;

        for(int k=a;k<b;k++){
            int p=example.simulated_particles(k);
            TV Vn_interpolate,V_pic,V_flip=particles.V(p);
            MATRIX<T,TV::m> B,grad_Vp;

            for(PARTICLE_GRID_ITERATOR<TV> it(example.weights,p,true,scratch);it.Valid();it.Next()){
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
                for(PARTICLE_GRID_ITERATOR<TV> it(example.weights,p,false,scratch);it.Valid();it.Next()){
                    TV_INT index=it.Index();
                    TV V_grid=example.velocity_new(index);
                    TV Z=example.grid.Center(index);
                    TV xi_new,xp_new;
                    if(example.use_midpoint){
                        xi_new=Z+dt/2*(V_grid+example.velocity(index));
                        xp_new=particles.X(p)+dt/2*(Vn_interpolate+V_pic);}
                    else{
                        xi_new=Z+dt*V_grid;
                        xp_new=particles.X(p)+dt*V_pic;}
                    B+=it.Weight()/2*(MATRIX<T,TV::m>::Outer_Product(V_grid,Z-particles.X(p)+xi_new-xp_new)
                            +MATRIX<T,TV::m>::Outer_Product(Z-particles.X(p)-xi_new+xp_new,V_grid));}

                    particles.V(p)=V_pic;
                    Perform_Particle_Collision(p);
                    if(example.use_midpoint) particles.X(p)+=(particles.V(p)+Vn_interpolate)*(dt/2);
                    else particles.X(p)+=particles.V(p)*dt;
                    particles.V(p)=V_flip*example.flip+V_pic*(1-example.flip);
                    particles.B(p)=B;

                    if(!example.grid.domain.Lazy_Inside(particles.X(p))) particles.valid(p)=false;}}
}
//#####################################################################
// Function Apply_Forces
//#####################################################################
template<class TV> void MPM_DRIVER<TV>::
Apply_Forces()
{
    example.Capture_Stress();
    example.Precompute_Forces(example.time);
    objective.Reset();
    if(example.use_forward_euler){
        objective.tmp2*=0;
        example.Add_Forces(objective.tmp2.u,example.time);
#pragma omp parallel for
        for(int i=0;i<example.valid_grid_indices.m;i++){
            int p=example.valid_grid_indices(i);
            dv.u.array(p)=example.dt/example.mass.array(p)*objective.tmp2.u.array(p);}}
    else{
        NEWTONS_METHOD<T> newtons_method;
        newtons_method.tolerance=example.newton_tolerance*example.dt;
        newtons_method.progress_tolerance=1e-5;
        newtons_method.max_iterations=example.newton_iterations;
        newtons_method.krylov_tolerance=example.solver_tolerance;
        newtons_method.max_krylov_iterations=example.solver_iterations;
        newtons_method.use_cg=true;
        newtons_method.debug=true;

        newtons_method.require_one_iteration=!objective.Initial_Guess(dv,newtons_method.tolerance);
        LOG::printf("max velocity: %P\n",Max_Particle_Speed());
        if(example.test_diff) objective.Test_Diff(dv);

        bool converged=newtons_method.Newtons_Method(objective,objective.system,dv,av);
        if(!converged) LOG::cout<<"WARNING: Newton's method did not converge"<<std::endl;}

    Apply_Friction();
    objective.Restore_F();

#pragma omp parallel for
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
    if(!example.use_particle_collision) return;;
    for(int i=0;i<example.collision_objects.m;i++){
        TV X=example.particles.X(p);
        IMPLICIT_OBJECT<TV>* io=example.collision_objects(i).io;
        T phi=io->Extended_Phi(X);
        if(phi>=0) continue;
        if(example.collision_objects(i).sticky) return;
        X-=phi*io->Normal(X);
        example.particles.X(p)=X;}
}
//#####################################################################
// Function Apply_Friction
//#####################################################################
template<class TV> void MPM_DRIVER<TV>::
Apply_Friction()
{
    if(!example.collision_objects.m) return;
    objective.Adjust_For_Collision(dv);
    objective.Compute_Unconstrained(dv,0,&objective.tmp0,0);
    objective.tmp1=objective.tmp0;
    objective.Project_Gradient_And_Prune_Constraints(objective.tmp1,true);

    objective.v1.u.array.Subset(objective.system.stuck_nodes).Fill(TV());
    for(int i=0;i<objective.system.collisions.m;i++){
        const typename MPM_KRYLOV_SYSTEM<TV>::COLLISION& c=objective.system.collisions(i);
        TV& v=objective.v1.u.array(c.p);
        T normal_force=TV::Dot_Product(c.n,objective.tmp0.u.array(c.p)-objective.tmp1.u.array(c.p));
        TV t=v.Projected_Orthogonal_To_Unit_Direction(c.n);
        T t_mag=t.Normalize();
        T coefficient_of_friction=example.collision_objects(c.object).friction;
        T k=coefficient_of_friction*normal_force/example.mass.array(c.p);
        if(t_mag<=k)
            v.Project_On_Unit_Direction(c.n);
        else v-=k*t;
        dv.u.array(c.p)=v-objective.v0.u.array(c.p);}
}
//#####################################################################
// Function Compute_Dt
//#####################################################################
template<class TV> typename TV::SCALAR MPM_DRIVER<TV>::
Compute_Dt() const
{
    T critical_speed=example.cfl*example.grid.DX().Min()/example.max_dt;
    T v=Grid_V_Upper_Bound();
    return (v>critical_speed)?(example.cfl*example.grid.DX().Min()/v):example.max_dt;
}
//#####################################################################
// Function Max_Particle_Speed
//#####################################################################
template<class TV> typename TV::SCALAR MPM_DRIVER<TV>::
Max_Particle_Speed() const
{
    T v2=0;
#pragma omp parallel for reduction(max:v2)
    for(int k=0;k<example.simulated_particles.m;k++){
        int p=example.simulated_particles(k);
        v2=max(v2,example.particles.V(p).Magnitude_Squared());}
    return sqrt(v2);
}
//#####################################################################
// Function Grid_V_Upper_Bound
//#####################################################################
template<class TV> typename TV::SCALAR MPM_DRIVER<TV>::
Grid_V_Upper_Bound() const
{
    if(!example.use_affine || !example.weights->constant_scalar_inertia_tensor) return Max_Particle_Speed();
    T result=0;
    T ksi=(T)6*sqrt((T)TV::m)*example.grid.One_Over_DX().Min();
#pragma omp parallel for reduction(max:result)
    for(int k=0;k<example.simulated_particles.m;k++){
        int p=example.simulated_particles(k);
        result=max(result,example.particles.V(p).Magnitude()+example.particles.B(p).Frobenius_Norm()*ksi);}
    return result;
}
//#####################################################################
// Function Update_Simulated_Particles
//#####################################################################
template<class TV> void MPM_DRIVER<TV>::
Update_Simulated_Particles()
{
    example.simulated_particles.Remove_All();
    for(int p=0;p<example.particles.number;p++)
        if(example.particles.valid(p))
            example.simulated_particles.Append(p);
    // example.simulated_particles=IDENTITY_ARRAY<>(example.particles.X.m);
    example.particle_is_simulated.Remove_All();
    example.particle_is_simulated.Resize(example.particles.X.m);
    example.particle_is_simulated.Subset(example.simulated_particles).Fill(true);

    for(int i=0;i<example.lagrangian_forces.m;i++){
        example.lagrangian_forces(i)->use_implicit_velocity_independent_forces=true;
        example.lagrangian_forces(i)->Update_Mpi(example.particle_is_simulated,0);}
}
//#####################################################################
// Function Print_Grid_Stats
//#####################################################################
template<class TV> void MPM_DRIVER<TV>::
Print_Grid_Stats(const char* str,T dt,const ARRAY<TV,TV_INT>& u,const ARRAY<TV,TV_INT>* u0)
{
    if(!example.print_stats) return;
    typename TV::SPIN am=example.Total_Grid_Angular_Momentum(dt,u,u0);
    TV lm=example.Total_Grid_Linear_Momentum(u);
    T ke=example.Total_Grid_Kinetic_Energy(u);
    LOG::cout<<str<<" linear  "<<"time " <<example.time<<" value "<<lm<<"  diff "<<(lm-example.last_linear_momentum)<<std::endl;
    LOG::cout<<str<<" angular "<<"time " <<example.time<<" value "<<am<<"  diff "<<(am-example.last_angular_momentum)<<std::endl;
    LOG::cout<<str<<" ke "<<"time " <<example.time<<" value "<<ke<<"  diff "<<(ke-example.last_grid_ke)<<std::endl;
    example.last_linear_momentum=lm;
    example.last_angular_momentum=am;
    example.last_grid_ke=ke;
}
//#####################################################################
// Function Print_Grid_Stats
//#####################################################################
template<class TV> void MPM_DRIVER<TV>::
Print_Particle_Stats(const char* str,T dt,const ARRAY<TV,TV_INT>& u,const ARRAY<TV,TV_INT>* u0)
{
    if(!example.print_stats) return;
    typename TV::SPIN am=example.Total_Particle_Angular_Momentum();
    TV lm=example.Total_Particle_Linear_Momentum();
    // T ke=example.Total_Particle_Kinetic_Energy();
    LOG::cout<<str<<" linear  "<<"time " <<example.time<<" value "<<lm<<"  diff "<<(lm-example.last_linear_momentum)<<std::endl;
    LOG::cout<<str<<" angular "<<"time " <<example.time<<" value "<<am<<"  diff "<<(am-example.last_angular_momentum)<<std::endl;
    // LOG::cout<<str<<" ke "<<ke<<"  diff "<<(ke-example.last_grid_ke)<<std::endl;
    example.last_linear_momentum=lm;
    example.last_angular_momentum=am;
    // example.last_grid_ke=ke;
}
//#####################################################################
// Function Print_Energy_Stats
//#####################################################################
template<class TV> void MPM_DRIVER<TV>::
Print_Energy_Stats(const char* str,const ARRAY<TV,TV_INT>& u)
{
    if(!example.print_stats) return;
    example.Capture_Stress();
    example.Precompute_Forces(example.time);
    T ke=example.Total_Grid_Kinetic_Energy(u);
    T pe=example.Potential_Energy(example.time);
    T te=ke+pe;
    LOG::cout<<str<<" kinetic  "<<"time " <<example.time<<" value "<<ke<<std::endl;
    LOG::cout<<str<<" potential "<<"time " <<example.time<<" value "<<pe<<std::endl;
    LOG::cout<<str<<" total energy "<<"time " <<example.time<<" value "<<te<<" diff "<<(te-example.last_te)<<std::endl;
    example.last_te=te;
}
//#####################################################################
namespace PhysBAM{
template class MPM_DRIVER<VECTOR<float,2> >;
template class MPM_DRIVER<VECTOR<float,3> >;
template class MPM_DRIVER<VECTOR<double,2> >;
template class MPM_DRIVER<VECTOR<double,3> >;
}
