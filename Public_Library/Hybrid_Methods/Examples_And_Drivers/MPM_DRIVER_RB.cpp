//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Log/DEBUG_SUBSTEPS.h>
#include <Core/Log/LOG.h>
#include <Core/Log/SCOPE.h>
#include <Tools/Krylov_Solvers/CONJUGATE_GRADIENT.h>
#include <Tools/Krylov_Solvers/MINRES.h>
#include <Tools/Nonlinear_Equations/NEWTONS_METHOD.h>
#include <Grid_Tools/Grids/CELL_ITERATOR.h>
#include <Grid_Tools/Grids/FACE_ITERATOR.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <Rigids/Collisions/COLLISION_HELPER.h>
#include <Deformables/Collisions_And_Interactions/IMPLICIT_OBJECT_COLLISION_PENALTY_FORCES.h>
#include <Deformables/Constitutive_Models/DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE.h>
#include <Solids/Solids/SOLID_BODY_COLLECTION.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_DRIVER_RB.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_EXAMPLE_RB.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_PARTICLES.h>
#include <Hybrid_Methods/Forces/MPM_FINITE_ELEMENTS.h>
#include <Hybrid_Methods/Forces/MPM_PLASTIC_FINITE_ELEMENTS.h>
#include <Hybrid_Methods/Forces/MPM_PLASTICITY_MODEL.h>
#include <Hybrid_Methods/Iterators/GATHER_SCATTER.h>
#include <Hybrid_Methods/Iterators/PARTICLE_GRID_WEIGHTS.h>
#include <Hybrid_Methods/System/MPM_KRYLOV_VECTOR.h>
#include <Hybrid_Methods/System/MPM_OBJECTIVE.h>
#ifdef USE_OPENMP
#include <omp.h>
#endif
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> MPM_DRIVER_RB<TV>::
MPM_DRIVER_RB(MPM_EXAMPLE_RB<TV>& example)
    :example(example),objective(*new MPM_OBJECTIVE<TV>(example)),
    dv(*new MPM_KRYLOV_VECTOR<TV>(example.valid_grid_indices)),
    rhs(*new MPM_KRYLOV_VECTOR<TV>(example.valid_grid_indices))
{
    DEBUG_SUBSTEPS::write_substeps_level=example.write_substeps_level;
    DEBUG_SUBSTEPS::writer=[=](const std::string& title){Write_Substep(title);};
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> MPM_DRIVER_RB<TV>::
~MPM_DRIVER_RB()
{
    DEBUG_SUBSTEPS::writer=0;
    delete &objective;
    delete &dv;
    delete &rhs;
    av.Delete_Pointers_And_Clean_Memory();
}
//#####################################################################
// Execute_Main_Program
//#####################################################################
template<class TV> void MPM_DRIVER_RB<TV>::
Execute_Main_Program()
{
    Initialize();
    Simulate_To_Frame(example.last_frame);
}
//#####################################################################
// Initialize
//#####################################################################
template<class TV> void MPM_DRIVER_RB<TV>::
Initialize()
{
    LOG::cout<<std::setprecision(16)<<std::endl;
    DEBUG_SUBSTEPS::write_substeps_level=example.substeps_delay_frame<0?example.write_substeps_level:-1;

    // setup time
    output_number=current_frame=example.restart;

    example.Initialize();
    PHYSBAM_ASSERT(example.grid.Is_MAC_Grid());
    if(example.restart)
        example.Read_Output_Files(example.restart);

    example.mass.Resize(example.grid.Domain_Indices(example.ghost));
    example.velocity.Resize(example.grid.Domain_Indices(example.ghost));
    example.velocity_save.Resize(example.grid.Domain_Indices(example.ghost));
    dv.u.Resize(example.grid.Domain_Indices(example.ghost));
    rhs.u.Resize(example.grid.Domain_Indices(example.ghost));
    objective.system.tmp.u.Resize(example.grid.Domain_Indices(example.ghost));

    example.particles.Store_B(example.use_affine);
    example.particles.Store_S(example.use_oldroyd);
    if(example.particles.store_B) example.Dp_inv.Resize(example.particles.X.m);

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

    // TODO
    // for(CELL_ITERATOR<TV> iterator(example.grid,example.ghost,GRID<TV>::GHOST_REGION);iterator.Valid();iterator.Next()){
    //     for(int k=0;k<fluid_walls.m;k++){
    //         if(fluid_walls(k)->Lazy_Inside(iterator.Location())){
    //             cell_soli....
    Update_Simulated_Particles();

    if(!example.restart) Write_Output_Files(0);
    PHYSBAM_DEBUG_WRITE_SUBSTEP("after init",1);
}
//#####################################################################
// Function Advance_One_Time_Step
//#####################################################################
template<class TV> void MPM_DRIVER_RB<TV>::
Advance_One_Time_Step()
{
    if(example.begin_time_step) example.begin_time_step(example.time);

    Update_Simulated_Particles();
    Print_Particle_Stats("particle state",example.dt);
    Update_Particle_Weights();
    example.gather_scatter.Prepare_Scatter(example.particles);
    Particle_To_Grid();
    Print_Grid_Stats("after particle to grid",example.dt,example.velocity,0);
    Print_Energy_Stats("after particle to grid",example.velocity);
    PHYSBAM_DEBUG_WRITE_SUBSTEP("after particle to grid",1);
    Apply_Forces();
    Print_Grid_Stats("after forces",example.dt,example.velocity,&example.velocity_save);
    PHYSBAM_DEBUG_WRITE_SUBSTEP("after forces",1);
    Grid_To_Particle_Limit_Dt();
    Limit_Dt_Sound_Speed();
    Grid_To_Particle();
    PHYSBAM_DEBUG_WRITE_SUBSTEP("after grid to particle",1);
    Update_Plasticity_And_Hardening();

    if(example.end_time_step) example.end_time_step(example.time);
}
//#####################################################################
// Function Simulate_To_Frame
//#####################################################################
template<class TV> void MPM_DRIVER_RB<TV>::
Simulate_To_Frame(const int frame)
{
    for(;current_frame<frame;current_frame++){
        LOG::SCOPE scope("FRAME","frame %d",current_frame+1);
        if(example.begin_frame) example.begin_frame(current_frame);
        if(example.substeps_delay_frame==current_frame)
            DEBUG_SUBSTEPS::write_substeps_level=example.write_substeps_level;
        T time_at_frame=example.time+example.frame_dt;
        bool done=false;
        for(int substep=0;!done;substep++){
            LOG::SCOPE scope("SUBSTEP","substep %d",substep+1);
            example.dt=Compute_Dt();
            example.dt=clamp(example.dt,example.min_dt,example.max_dt);
            T next_time=example.time+example.dt;
            if(next_time+(T).0001*example.dt>time_at_frame){ // Account for roundoff errors
                next_time=time_at_frame;
                done=true;}
            else if(next_time+example.dt>time_at_frame) next_time=(example.time+time_at_frame)/2;
            example.dt=next_time-example.time;
            LOG::cout<<"substep dt: "<<example.dt<<std::endl;

            Advance_One_Time_Step();

            // Time step was reduced
            if(example.dt<next_time-example.time){
                LOG::printf("dt reduced: %g to %g\n",next_time-example.time,example.dt);
                next_time=example.time+example.dt;
                done=false;}

            PHYSBAM_DEBUG_WRITE_SUBSTEP("end substep %i",0,substep);
            example.time=next_time;}
        if(example.end_frame) example.end_frame(current_frame);
        Write_Output_Files(++output_number);}
}
//#####################################################################
// Function Write_Substep
//#####################################################################
template<class TV> void MPM_DRIVER_RB<TV>::
Write_Substep(const std::string& title)
{
    example.frame_title=title;
    LOG::printf("Writing substep [%s]: output_number=%i, time=%g, frame=%i\n",
        example.frame_title,output_number+1,example.time,current_frame);
    Write_Output_Files(++output_number);
    example.frame_title="";
}
//#####################################################################
// Write_Output_Files
//#####################################################################
template<class TV> void MPM_DRIVER_RB<TV>::
Write_Output_Files(const int frame)
{
    LOG::SCOPE scope("Write_Output_Files");
    Create_Directory(example.output_directory);
    Create_Directory(example.output_directory+LOG::sprintf("/%d",frame));
    Create_Directory(example.output_directory+"/common");
    Write_To_Text_File(example.output_directory+LOG::sprintf("/%d/frame_title",frame),example.frame_title);
    if(frame==0)
        Write_To_Text_File(example.output_directory+"/common/first_frame",frame,"\n");
    example.Write_Output_Files(frame);
    Write_To_Text_File(example.output_directory+"/common/last_frame",frame,"\n");
}
//#####################################################################
// Function Update_Particle_Weights
//#####################################################################
template<class TV> void MPM_DRIVER_RB<TV>::
Update_Particle_Weights()
{
    example.weights->Update(example.particles.X);
    if(example.particles.store_B) example.weights->Dp_Inverse(example.particles.X,example.Dp_inv);
}
//#####################################################################
// Function Particle_To_Grid
//#####################################################################
template<class TV> void MPM_DRIVER_RB<TV>::
Particle_To_Grid()
{
    MPM_PARTICLES<TV>& particles=example.particles;

#pragma omp parallel for
    for(int i=0;i<example.mass.array.m;i++){
        example.mass.array(i)=0;
        example.velocity.array(i)=TV();
        example.velocity_save.array(i)=TV();}

    bool use_gradient=example.weights->use_gradient_transfer;
    example.gather_scatter.template Scatter<int>(true,
        [this,&particles,use_gradient](int p,const PARTICLE_GRID_ITERATOR<TV>& it,int data)
        {
            T w=it.Weight();
            TV_INT index=it.Index();
            example.mass(index)+=w*particles.mass(p);
            TV V=w*particles.V(p);
            if(example.use_affine){
                TV Y=example.grid.Center(index)-particles.X(p);
                if(!example.lag_Dp) Y=example.Dp_inv(p)*Y;
                V+=particles.B(p)*Y*w;}
            example.velocity(index)+=particles.mass(p)*V;
        });

    example.valid_grid_indices.Remove_All();
    example.valid_grid_cell_indices.Remove_All();

    Reflection_Boundary_Condition(example.velocity,true);
    Reflection_Boundary_Condition(example.mass,false);

    for(RANGE_ITERATOR<TV::m> it(example.mass.domain);it.Valid();it.Next()){
        int i=example.mass.Standard_Index(it.index);
        if(example.mass.array(i)){
            example.valid_grid_indices.Append(i);
            example.valid_grid_cell_indices.Append(it.index);
            TV v=example.velocity.array(i)/example.mass.array(i);
            example.velocity.array(i)=v;
            example.velocity_save.array(i)=v;}
        else example.velocity.array(i)=TV();}

    Rasterize_Rigid_Bodies();
}
//#####################################################################
// Function Grid_To_Particle
//#####################################################################
template<class TV> void MPM_DRIVER_RB<TV>::
Grid_To_Particle()
{
    struct HELPER
    {
        TV V_pic,V_weight_old,V_pic_fric;
        MATRIX<T,TV::m> B,grad_Vp,grad_Vp_fric;
    };

    MPM_PARTICLES<TV>& particles=example.particles;
    T dt=example.dt;
    bool use_gradient_transfer=example.weights->use_gradient_transfer;
    bool need_dv_fric=(example.use_affine && use_gradient_transfer);

    example.gather_scatter.template Gather<HELPER>(true,
        [](int p,HELPER& h){h=HELPER();},
        [this,need_dv_fric,dt,&particles](int p,const PARTICLE_GRID_ITERATOR<TV>& it,HELPER& h)
        {
            T w=it.Weight();
            TV dw=it.Gradient();
            TV_INT index=it.Index();
            TV V_new=example.velocity(index);
            TV V_old=example.velocity_save(index);
            TV V_fric=example.velocity_friction_save(index);
            h.V_pic+=w*V_new;
            h.V_weight_old+=w*V_old;
            h.V_pic_fric+=w*V_fric;
            TV V_grid=example.use_midpoint?(T).5*(V_new+V_old):V_new;
            h.grad_Vp+=Outer_Product(V_grid,dw);
            if(need_dv_fric){
                TV V_grid_fric=example.use_midpoint?(T).5*(V_fric+V_old):V_fric;
                h.grad_Vp_fric+=Outer_Product(V_grid_fric,dw);}
            if(example.use_affine && !example.lag_Dp){
                TV Z=example.grid.Center(index);
                TV xi_new=Z+dt*V_grid;
                h.B+=w/2*(Outer_Product(V_fric,Z+xi_new)+Outer_Product(Z-xi_new,V_fric));}
        },
        [this,dt,&particles](int p,HELPER& h)
        {
            MATRIX<T,TV::m> A=dt*h.grad_Vp+1;
            if(example.quad_F_coeff) A+=sqr(dt*h.grad_Vp)*example.quad_F_coeff;
            particles.F(p)=A*particles.F(p);
            TV xp_new;
            if(example.use_midpoint) xp_new=particles.X(p)+dt/2*(h.V_weight_old+h.V_pic);
            else xp_new=particles.X(p)+dt*h.V_pic;
            if(particles.store_S){
                T k=example.dt*example.inv_Wi;
                particles.S(p)=(SYMMETRIC_MATRIX<T,TV::m>::Conjugate(A,particles.S(p))+k)/(1+k);}
            if(example.use_affine){
                h.B-=(T).5*(Outer_Product(h.V_pic_fric,particles.X(p)+xp_new)
                    +Outer_Product(particles.X(p)-xp_new,h.V_pic_fric));
                if(example.weights->use_gradient_transfer) particles.B(p)=h.grad_Vp_fric;
                else if(example.lag_Dp) particles.B(p)=h.B*example.Dp_inv(p);
                else particles.B(p)=h.B;}
            particles.X(p)=xp_new;
            TV V_flip=particles.V(p)+h.V_pic_fric-h.V_weight_old;
            particles.V(p)=V_flip*example.flip+h.V_pic_fric*(1-example.flip);

            Reflect_Or_Invalidate_Particle(p);
        });

    Move_Rigid_Bodies();
}
// Lower s if necessary so that |t*(a+t*b)|<=bound for all 0<=t<=s
template<class T>
void Enforce_Limit_Max(T& s,T bound,T a,T b)
{
    if(a<0){a=-a;b=-b;}
    if(a*a+4*b*bound>=0){
        if(s*(a+s*b)>bound)
            s=(2*bound)/(a+sqrt(a*a+4*b*bound));}
    else{
        if(s*(a+s*b)<-bound)
            s=-(a+sqrt(a*a-4*b*bound))/(2*b);}
}
// Lower s if necessary so that t*(a+t*b)<=bound for all 0<=t<=s
template<class T,int d>
void Enforce_Limit_Max(T& s,T bound,const VECTOR<T,d>& a,const VECTOR<T,d>& b)
{
    for(int i=0;i<d;i++) Enforce_Limit_Max(s,bound,a(i),b(i));
}
// Lower s if necessary so that t*(a+t*b)<=bound for all 0<=t<=s
template<class T,int d>
void Enforce_Limit_Max(T& s,T bound,const MATRIX<T,d>& a,const MATRIX<T,d>& b)
{
    for(int i=0;i<d;i++)
        for(int j=0;j<d;j++)
            Enforce_Limit_Max(s,bound,a(i,j),b(i,j));
}
//#####################################################################
// Function Grid_To_Particle_Limit_Dt
//#####################################################################
template<class TV> void MPM_DRIVER_RB<TV>::
Grid_To_Particle_Limit_Dt()
{
    if(!example.use_strong_cfl) return;
    struct HELPER
    {
        TV V_pic,V_pic_s,V_weight_old;
        MATRIX<T,TV::m> grad_Vp,grad_Vp_s;
    };

    T dt=example.dt,s=1;
    T midpoint_frac=example.use_midpoint?(T).5:1;

    // TODO: this is NOT threadsafe.
    example.gather_scatter.template Gather<HELPER>(true,
        [](int p,HELPER& h){h=HELPER();},
        [this,midpoint_frac](int p,const PARTICLE_GRID_ITERATOR<TV>& it,HELPER& h)
        {
            T w=it.Weight();
            TV dw=it.Gradient();
            TV_INT index=it.Index();
            TV V_new=example.velocity(index);
            TV V_old=example.velocity_save(index);
            TV V_fric=example.velocity_friction_save(index);

            h.V_pic+=w*V_old;
            h.V_pic_s+=w*(V_new-V_old);
            h.V_weight_old+=w*V_old;
            TV V_grid=V_old;
            TV V_grid_s=midpoint_frac*(V_new-V_old);
            h.grad_Vp+=Outer_Product(V_grid,dw);
            h.grad_Vp_s+=Outer_Product(V_grid_s,dw);
        },
        [this,dt,&s](int p,HELPER& h)
        {
            Enforce_Limit_Max(s,example.cfl_F,dt*h.grad_Vp,dt*h.grad_Vp_s);
            TV xp_new_s,xp_new_s2;
            if(example.use_midpoint){xp_new_s=dt/2*(h.V_weight_old+h.V_pic);xp_new_s2=dt/2*h.V_pic_s;}
            else{xp_new_s=dt*h.V_pic;xp_new_s2=dt*h.V_pic_s;}
            Enforce_Limit_Max(s,example.cfl,xp_new_s,xp_new_s2);
        });
    if(example.dt*s<example.min_dt) s=example.min_dt/example.dt;
    if(s>=1) return;
    LOG::printf("X J CFL scale: %g -> %g\n",example.dt,example.dt*s);
    example.dt*=s;
    example.velocity.array=(example.velocity.array-example.velocity_save.array)*s+example.velocity_save.array;
    example.velocity_friction_save.array=(example.velocity_friction_save.array-example.velocity_save.array)*s+example.velocity_save.array;
}
//#####################################################################
// Function Limit_Dt_Sound_Speed
//#####################################################################
template<class TV> void MPM_DRIVER_RB<TV>::
Limit_Dt_Sound_Speed()
{
    if(!example.use_sound_speed_cfl) return;
    T dt=example.dt;
    for(int f=0;f<example.forces.m;f++){
        if(const MPM_FINITE_ELEMENTS<TV>* force=dynamic_cast<MPM_FINITE_ELEMENTS<TV>*>(example.forces(f))){
            for(int k=0;k<example.simulated_particles.m;k++){
                int p=example.simulated_particles(k);
                const DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<TV>& d=force->dPi_dF(p);
                T K=d.H.Diagonal_Part().x.Max_Abs();
                T J=example.particles.F(p).Determinant();
                T rho=example.particles.mass(p)/(example.particles.volume(p)*J);
                T speed=sqrt(K/(J*rho));
                T new_dt=example.grid.dX.Min()/speed*example.cfl_sound;
                dt=min(dt,new_dt);}}}
    if(dt<example.min_dt) dt=example.min_dt;
    if(dt>=example.dt) return;
    T s=dt/example.dt;
    LOG::printf("SOUND CFL %g %g (%g)\n",example.dt,dt,s);
    example.velocity.array=(example.velocity.array-example.velocity_save.array)*s+example.velocity_save.array;
    example.velocity_friction_save.array=(example.velocity_friction_save.array-example.velocity_save.array)*s+example.velocity_save.array;
    example.dt=dt;
}
//#####################################################################
// Function Update_Plasticity_And_Hardening
//#####################################################################
template<class TV> void MPM_DRIVER_RB<TV>::
Update_Plasticity_And_Hardening()
{
    for(int i=0;i<example.plasticity_models.m;i++)
        example.plasticity_models(i)->Update_Particles();
}
//#####################################################################
// Function Apply_Forces
//#####################################################################
template<class TV> void MPM_DRIVER_RB<TV>::
Apply_Forces()
{
    example.Capture_Stress();
    objective.Reset();
    LOG::printf("max velocity: %.16P\n",Max_Particle_Speed());
    if(example.use_symplectic_euler){
        example.Precompute_Forces(example.time,example.dt,0);
        objective.tmp2.u*=0;
        example.Add_Forces(objective.tmp2.u,example.time);
        Reflection_Boundary_Condition(objective.tmp2.u,true);
        Apply_Rigid_Body_Forces();
        if(example.pairwise_collisions){
#pragma omp parallel for
            for(int i=0;i<example.valid_grid_indices.m;i++){
                int p=example.valid_grid_indices(i);
                example.velocity.array(p)+=example.dt/example.mass.array(p)*objective.tmp2.u.array(p);}
            Process_Pairwise_Collisions();
            example.velocity_friction_save=example.velocity;}
        else{
            for(int i=0;i<example.valid_grid_indices.m;i++){
                int p=example.valid_grid_indices(i);
                dv.u.array(p)=example.dt/example.mass.array(p)*objective.tmp2.u.array(p);}
            objective.system.forced_collisions.Remove_All();
            objective.tmp1=dv;
            objective.Adjust_For_Collision(dv);
            objective.tmp0=dv;
            Apply_Friction();}}
    else{
        NEWTONS_METHOD<T> newtons_method;
        newtons_method.tolerance=example.newton_tolerance*example.dt;
        newtons_method.progress_tolerance=1e-5;
        newtons_method.max_iterations=example.newton_iterations;
        newtons_method.krylov_tolerance=example.solver_tolerance;
        newtons_method.max_krylov_iterations=example.solver_iterations;
        newtons_method.use_cg=true;
        newtons_method.debug=true;
        if(example.asymmetric_system){
            newtons_method.use_gmres=true;
            newtons_method.use_cg=false;
            newtons_method.Make_Vanilla_Newton();}
        
        example.Update_Lagged_Forces(example.time);
        newtons_method.require_one_iteration=!objective.Initial_Guess(dv,newtons_method.tolerance,example.asymmetric_system);
        if(example.test_diff) objective.Test_Diff(dv);

        objective.system.forced_collisions.Remove_All();

        bool converged=newtons_method.Newtons_Method(objective,objective.system,dv,av);
        Apply_Rigid_Body_Forces();

        if(!converged) LOG::cout<<"WARNING: Newton's method did not converge"<<std::endl;
        Apply_Friction();
        objective.Restore_F();}

    if(!example.pairwise_collisions){
#pragma omp parallel for
        for(int i=0;i<example.valid_grid_indices.m;i++){
            int j=example.valid_grid_indices(i);
            example.velocity.array(j)=dv.u.array(j)+objective.v0.u.array(j);
            example.velocity_friction_save.array(j)+=objective.v0.u.array(j);}
        example.velocity_friction_save.array.Subset(objective.system.stuck_nodes)=objective.system.stuck_velocity;}
}
//#####################################################################
// Function Apply_Friction
//#####################################################################
template<class TV> void MPM_DRIVER_RB<TV>::
Apply_Friction()
{
    example.velocity_friction_save=dv.u;
    if(!example.collision_objects.m) return;
    if(example.use_symplectic_euler){
        objective.v1.Copy(1,objective.v0,dv);}
    else{
        objective.Compute_Unconstrained(dv,0,&objective.tmp0,0);
        objective.tmp1=objective.tmp0;
        objective.Project_Gradient_And_Prune_Constraints(objective.tmp1,true);}

    objective.v1.u.array.Subset(objective.system.stuck_nodes).Fill(TV());
    for(int i=0;i<objective.system.collisions.m;i++){
        const typename MPM_KRYLOV_SYSTEM<TV>::COLLISION& c=objective.system.collisions(i);
        TV& v=objective.v1.u.array(c.p);
        // Note: tmp0,tmp1 already have mass removed; they have units of velocity
        T normal_force=TV::Dot_Product(c.n,objective.tmp0.u.array(c.p)-objective.tmp1.u.array(c.p));
        TV t=v.Projected_Orthogonal_To_Unit_Direction(c.n);
        T t_mag=t.Normalize();
        T coefficient_of_friction=example.collision_objects(c.object)->friction;
        T k=coefficient_of_friction*normal_force;
        if(t_mag<=k)
            v.Project_On_Unit_Direction(c.n);
        else v-=k*t;
        example.velocity_friction_save.array(c.p)=v-objective.v0.u.array(c.p);}
}
//#####################################################################
// Function Compute_Dt
//#####################################################################
template<class TV> typename TV::SCALAR MPM_DRIVER_RB<TV>::
Compute_Dt() const
{
    T critical_speed=example.cfl*example.grid.dX.Min()/example.max_dt;
    T v=Grid_V_Upper_Bound();
    return (v>critical_speed)?(example.cfl*example.grid.dX.Min()/v):example.max_dt;
}
//#####################################################################
// Function Max_Particle_Speed
//#####################################################################
template<class TV> typename TV::SCALAR MPM_DRIVER_RB<TV>::
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
template<class TV> typename TV::SCALAR MPM_DRIVER_RB<TV>::
Grid_V_Upper_Bound() const
{
    if(!example.use_affine) return Max_Particle_Speed();
    T result=0;
    T xi=(T)6*sqrt((T)TV::m)*example.grid.one_over_dX.Min();
#pragma omp parallel for reduction(max:result)
    for(int k=0;k<example.simulated_particles.m;k++){
        int p=example.simulated_particles(k);
        T v=example.particles.V(p).Magnitude();
        if(example.particles.store_B) v+=example.particles.B(p).Frobenius_Norm()*xi;
        result=max(result,v);}
    return result;
}
//#####################################################################
// Function Update_Simulated_Particles
//#####################################################################
template<class TV> void MPM_DRIVER_RB<TV>::
Update_Simulated_Particles()
{
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=example.solid_body_collection.rigid_body_collection;
    RIGID_BODY_PARTICLES<TV>& rigid_body_particles=rigid_body_collection.rigid_body_particles;

    example.simulated_particles.Remove_All();
    for(int p=0;p<example.particles.number;p++)
        if(example.particles.valid(p))
            example.simulated_particles.Append(p);
    // example.simulated_particles=IDENTITY_ARRAY<>(example.particles.X.m);
    example.particle_is_simulated.Remove_All();
    example.particle_is_simulated.Resize(example.particles.X.m);
    example.particle_is_simulated.Subset(example.simulated_particles).Fill(true);

    for(int i=0;i<example.lagrangian_forces.m;i++)
        example.lagrangian_forces(i)->Update_Mpi(example.particle_is_simulated,0);

    example.rigid_body_is_simulated.Resize(rigid_body_particles.frame.m);
    for(int b=0;b<rigid_body_particles.frame.m;b++)
        example.rigid_body_is_simulated(b)=rigid_body_collection.Is_Active(b);

    for(int i=0;i<rigid_body_collection.rigids_forces.m;i++)
        rigid_body_collection.rigids_forces(i)->Update_Mpi(example.rigid_body_is_simulated);
}
//#####################################################################
// Function Print_Grid_Stats
//#####################################################################
template<class TV> void MPM_DRIVER_RB<TV>::
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
template<class TV> void MPM_DRIVER_RB<TV>::
Print_Particle_Stats(const char* str,T dt)
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
template<class TV> void MPM_DRIVER_RB<TV>::
Print_Energy_Stats(const char* str,const ARRAY<TV,TV_INT>& u)
{
    if(!example.print_stats) return;
    example.Capture_Stress();
    example.Precompute_Forces(example.time,example.dt,false);
    T ke=example.Total_Grid_Kinetic_Energy(u);
    T ke2=example.Total_Particle_Kinetic_Energy();
    T pe=example.Potential_Energy(example.time);
    T te=ke+pe;
    LOG::cout<<str<<" kinetic  "<<"time " <<example.time<<" value "<<ke<<std::endl;
    LOG::cout<<str<<" potential "<<"time " <<example.time<<" value "<<pe<<std::endl;
    LOG::cout<<str<<" total energy "<<"time " <<example.time<<" value "<<te<<" diff "<<(te-example.last_te)<<std::endl;
    LOG::cout<<str<<" particle total energy "<<"time " <<example.time<<" value "<<(ke2+pe)<<std::endl;
    example.last_te=te;
}
template<class T,int d> void flip(VECTOR<T,d>& u,int a){u(a)=-u(a);}
template<class T> void flip(T& u,int a){}
//#####################################################################
// Function Reflection_Boundary_Condition
//#####################################################################
template<class TV> template<class S> void MPM_DRIVER_RB<TV>::
Reflection_Boundary_Condition(ARRAY<S,TV_INT>& u,bool flip_sign)
{
    if(!example.reflection_bc_flags) return;
    TV_INT ranges[2]={TV_INT(),example.grid.numbers_of_cells};
    for(RANGE_ITERATOR<TV::m> it(example.grid.Domain_Indices(),example.ghost,0,
            RI::ghost|RI::side_mask,example.reflection_bc_flags);it.Valid();it.Next()){
        int axis=it.side/2,side=it.side%2;
        TV_INT pair=it.index;
        pair(axis)=2*ranges[side](axis)-it.index(axis)-1;
        if(flip_sign) flip(u(it.index),axis);
        u(pair)+=u(it.index);
        u(it.index)=S();}
    for(RANGE_ITERATOR<TV::m> it(example.grid.Domain_Indices(),example.ghost,0,
            RI::ghost|RI::delay_corners|RI::side_mask,example.reflection_bc_flags);it.Valid();it.Next()){
        int axis=it.side/2,side=it.side%2;
        TV_INT pair=it.index;
        pair(axis)=2*ranges[side](axis)-it.index(axis)-1;
        u(it.index)=u(pair);
        if(flip_sign) flip(u(it.index),axis);}
}
//#####################################################################
// Function Reflect_Particles
//#####################################################################
template<class TV> void MPM_DRIVER_RB<TV>::
Reflect_Or_Invalidate_Particle(int p)
{
    int f=example.reflection_bc_flags;
    TV A=example.grid.domain.min_corner;
    TV B=example.grid.domain.max_corner;
    TV& X=example.particles.X(p),&V=example.particles.V(p);
    for(int a=0;a<TV::m;a++){
        if(X(a)<A(a)){
            if((f>>(2*a))&1){X(a)=2*A(a)-X(a);V(a)=-V(a);}
            else example.particles.valid(p)=false;}
        else if(X(a)>B(a)){
            if((f>>(2*a+1))&1){X(a)=2*B(a)-X(a);V(a)=-V(a);}
            else example.particles.valid(p)=false;}}
}
//#####################################################################
// Function Update_Rotation_Helper
//#####################################################################
template<class T>
inline VECTOR<T,3> Update_Rotation_Helper(const T dt,RIGID_BODY<VECTOR<T,3> >& rigid_body)
{
    const VECTOR<T,3> &w=rigid_body.Twist().angular,&L=rigid_body.Angular_Momentum();
    return w-(T).5*dt*rigid_body.World_Space_Inertia_Tensor_Inverse_Times(w.Cross(L));
}
//#####################################################################
// Function Update_Rotation_Helper
//#####################################################################
template<class T>
inline VECTOR<T,1> Update_Rotation_Helper(const T dt,RIGID_BODY<VECTOR<T,2> >& rigid_body)
{
    return rigid_body.Twist().angular;
}
//#####################################################################
// Function Move_Rigid_Bodies
//#####################################################################
template<class TV> void MPM_DRIVER_RB<TV>::
Move_Rigid_Bodies()
{
    RIGID_BODY_PARTICLES<TV>& rigid_body_particles=example.solid_body_collection.rigid_body_collection.rigid_body_particles;
    for(int b=0;b<rigid_body_particles.frame.m;b++){
        if(!example.solid_body_collection.rigid_body_collection.Is_Active(b)) continue;
        rigid_body_particles.frame(b).t+=example.dt*rigid_body_particles.twist(b).linear;
        ROTATION<TV>& R=rigid_body_particles.frame(b).r;
        typename TV::SPIN rotate_amount=Update_Rotation_Helper(example.dt,*rigid_body_particles.rigid_body(b));
        R=ROTATION<TV>::From_Rotation_Vector(example.dt*rotate_amount)*R;
        R.Normalize();
        rigid_body_particles.rigid_body(b)->Update_Angular_Velocity();} // Note that the value of w changes here.
}
//#####################################################################
// Function Apply_Rigid_Body_Forces
//#####################################################################
template<class TV> void MPM_DRIVER_RB<TV>::
Apply_Rigid_Body_Forces()
{
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=example.solid_body_collection.rigid_body_collection;
    RIGID_BODY_PARTICLES<TV>& rigid_body_particles=rigid_body_collection.rigid_body_particles;
    rigid_forces.Resize(rigid_body_particles.frame.m);
    rigid_forces.Fill(TWIST<TV>());

    rigid_body_collection.Update_Position_Based_State(example.time+example.dt);
    rigid_body_collection.Add_Velocity_Independent_Forces(rigid_forces,example.time);
    
    for(int b=0;b<rigid_body_particles.frame.m;b++){
        if(!rigid_body_collection.Is_Active(b)) continue;
        rigid_body_particles.twist(b)+=rigid_body_particles.rigid_body(b)->Inertia_Inverse_Times(rigid_forces(b)*example.dt);}

    rigid_body_collection.Update_Angular_Momentum();
}
//#####################################################################
// Function Rasterize_Rigid_Bodies
//#####################################################################
template<class TV> void MPM_DRIVER_RB<TV>::
Rasterize_Rigid_Bodies()
{
    RIGID_BODY_PARTICLES<TV>& rigid_body_particles=example.solid_body_collection.rigid_body_collection.rigid_body_particles;
    T padding=example.grid.dX.Max()*2;
    RANGE<TV_INT> cell_domain=example.grid.Domain_Indices(example.ghost);
    example.rasterized_data.Resize(cell_domain);
    for(int b=0;b<rigid_body_particles.frame.m;b++){
        RIGID_BODY<TV>& rigid_body=*rigid_body_particles.rigid_body(b);
        rigid_body.Update_Bounding_Box_From_Implicit_Geometry();
        RANGE<TV> box=rigid_body.Axis_Aligned_Bounding_Box().Thickened(padding);
        RANGE<TV_INT> grid_range=example.grid.Clamp_To_Cell(box,example.ghost+1).Intersect(cell_domain);
        for(RANGE_ITERATOR<TV::m> it(grid_range);it.Valid();it.Next()){
            TV X=example.grid.Center(it.index);
            T phi=rigid_body.Implicit_Geometry_Extended_Value(X);
            if(phi<padding)
                example.rasterized_data.Insert(it.index,{b,phi});}}
}
//#####################################################################
// Function Process_Pairwise_Collisions
//#####################################################################
template<class TV> void MPM_DRIVER_RB<TV>::
Process_Pairwise_Collisions()
{
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=example.solid_body_collection.rigid_body_collection;
    bool need_another_iteration=true;
    for(int i=0;i<example.collision_iterations && need_another_iteration;i++){
        need_another_iteration=false;
        for(CELL_ITERATOR<TV> it(example.grid);it.Valid();it.Next()){
            TV X=it.Location();
            auto collisions=example.rasterized_data.Get(it.index);
            for(int c1=0;c1<collisions.m-1;c1++)
                for(int c2=c1+1;c2<collisions.m;c2++){
                    int r1=collisions(c1).id,r2=collisions(c2).id;
                    if(collisions(c1).phi+collisions(c2).phi<=0){
                        RIGID_BODY<TV>& body1=rigid_body_collection.Rigid_Body(r1);
                        RIGID_BODY<TV>& body2=rigid_body_collection.Rigid_Body(r2);
                        if(body1.Has_Infinite_Inertia() && body2.Has_Infinite_Inertia())
                            continue;
                        TV normal=(body2.Implicit_Geometry_Normal(X)-body1.Implicit_Geometry_Normal(X)).Normalized();
                        TV linear_impulse=PhysBAM::Compute_Collision_Impulse(
                            normal,
                            RIGID_BODY<TV>::Impulse_Factor(body1,body2,X),
                            RIGID_BODY<TV>::Relative_Velocity(body1,body2,X),
                            RIGID_BODY<TV>::Coefficient_Of_Restitution(body1,body2),
                            RIGID_BODY<TV>::Coefficient_Of_Friction(body1,body2),
                            0);
                        RIGID_BODY<TV>::Apply_Impulse(body1,body2,X,linear_impulse);
                        need_another_iteration=true;}}
            if(example.mass(it.index)){
                for(int c=0;c<collisions.m;c++){
                    if(collisions(c).phi>0) continue;
                    RIGID_BODY<TV>& body=rigid_body_collection.Rigid_Body(collisions(c).id);
                    TV linear_impulse=PhysBAM::Compute_Collision_Impulse(
                        -body.Implicit_Geometry_Normal(X),
                        body.Impulse_Factor(X)+1/example.mass(it.index),
                        body.Pointwise_Object_Velocity(X)-example.velocity(it.index),
                        body.coefficient_of_restitution,
                        body.coefficient_of_friction,
                        0);
                    body.Apply_Impulse_To_Body(X,linear_impulse);
                    example.velocity(it.index)-=linear_impulse/example.mass(it.index);
                    need_another_iteration=true;}}}

        // TODO: Do we actually need this?
        // force it to the last iteration (so it picks up contact pairs for rigid/rigid one time at least)
        if(!need_another_iteration && i<example.collision_iterations-1){
            i=example.collision_iterations-2;
            need_another_iteration=true;}}
}
//#####################################################################
namespace PhysBAM{
template class MPM_DRIVER_RB<VECTOR<float,2> >;
template class MPM_DRIVER_RB<VECTOR<float,3> >;
template class MPM_DRIVER_RB<VECTOR<double,2> >;
template class MPM_DRIVER_RB<VECTOR<double,3> >;
}
