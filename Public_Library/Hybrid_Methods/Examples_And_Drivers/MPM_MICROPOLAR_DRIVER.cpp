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
#include <Deformables/Collisions_And_Interactions/IMPLICIT_OBJECT_COLLISION_PENALTY_FORCES.h>
#include <Deformables/Constitutive_Models/DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_MICROPOLAR_DRIVER.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_MICROPOLAR_EXAMPLE.h>
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
template<class TV> MPM_MICROPOLAR_DRIVER<TV>::
MPM_MICROPOLAR_DRIVER(MPM_MICROPOLAR_EXAMPLE<TV>& example)
    :example(example)
{
    DEBUG_SUBSTEPS::write_substeps_level=example.write_substeps_level;
    DEBUG_SUBSTEPS::writer=[=](const std::string& title){Write_Substep(title);};
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> MPM_MICROPOLAR_DRIVER<TV>::
~MPM_MICROPOLAR_DRIVER()
{
    DEBUG_SUBSTEPS::writer=0;
}
//#####################################################################
// Execute_Main_Program
//#####################################################################
template<class TV> void MPM_MICROPOLAR_DRIVER<TV>::
Execute_Main_Program()
{
    Initialize();
    Simulate_To_Frame(example.last_frame);
}
//#####################################################################
// Initialize
//#####################################################################
template<class TV> void MPM_MICROPOLAR_DRIVER<TV>::
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
    forces.Resize(example.grid.Domain_Indices(example.ghost));

    example.particles.Store_C(example.use_affine);

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

    Update_Simulated_Particles();

    if(!example.restart) Write_Output_Files(0);
    PHYSBAM_DEBUG_WRITE_SUBSTEP("after init",1);
}
//#####################################################################
// Function Advance_One_Time_Step
//#####################################################################
template<class TV> void MPM_MICROPOLAR_DRIVER<TV>::
Advance_One_Time_Step()
{
    if(example.begin_time_step) example.begin_time_step(example.time);

    Update_Simulated_Particles();
    Update_Particle_Weights();
    example.gather_scatter.Prepare_Scatter(example.particles);
    Particle_To_Grid();
    PHYSBAM_DEBUG_WRITE_SUBSTEP("after particle to grid",1);
    Apply_Forces();
    PHYSBAM_DEBUG_WRITE_SUBSTEP("after forces",1);
    Grid_To_Particle_Limit_Dt();
    Limit_Dt_Sound_Speed();
    PHYSBAM_DEBUG_WRITE_SUBSTEP("after limits",1);
    Grid_To_Particle();
    PHYSBAM_DEBUG_WRITE_SUBSTEP("after grid to particle",1);

    if(example.end_time_step) example.end_time_step(example.time);
}
//#####################################################################
// Function Simulate_To_Frame
//#####################################################################
template<class TV> void MPM_MICROPOLAR_DRIVER<TV>::
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
template<class TV> void MPM_MICROPOLAR_DRIVER<TV>::
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
template<class TV> void MPM_MICROPOLAR_DRIVER<TV>::
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
template<class TV> void MPM_MICROPOLAR_DRIVER<TV>::
Update_Particle_Weights()
{
    example.weights->Update(example.particles.X);
}
//#####################################################################
// Function Particle_To_Grid
//#####################################################################
template<class TV> void MPM_MICROPOLAR_DRIVER<TV>::
Particle_To_Grid()
{
    MPM_PARTICLES<TV>& particles=example.particles;

#pragma omp parallel for
    for(int i=0;i<example.mass.array.m;i++){
        example.mass.array(i)=0;
        example.velocity.array(i)=TV();}

    example.gather_scatter.template Scatter<int>(true,
        [this,&particles](int p,const PARTICLE_GRID_ITERATOR<TV>& it,int data)
        {
            T w=it.Weight();
            TV_INT index=it.Index();
            example.mass(index)+=w*particles.mass(p);
            TV V=w*particles.V(p);
            if(example.use_affine){
                TV Y=example.grid.Center(index)-particles.X(p);
                V+=particles.C(p)*Y*w;}
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
            example.velocity.array(i)/=example.mass.array(i);}
        else example.velocity.array(i)=TV();}

    example.velocity_save=example.velocity;
}
//#####################################################################
// Function Grid_To_Particle
//#####################################################################
template<class TV> void MPM_MICROPOLAR_DRIVER<TV>::
Grid_To_Particle()
{
    struct HELPER
    {
        TV V_pic,V_weight_old;
        MATRIX<T,TV::m> B,grad_Vp;
    };

    MPM_PARTICLES<TV>& particles=example.particles;
    T dt=example.dt;
    ARRAY<SYMMETRIC_MATRIX<T,TV::m> > Dp_inv;
    if(example.use_affine && !example.weights->use_gradient_transfer){
        Dp_inv.Resize(example.particles.X.m);
        example.weights->Dp_Inverse(example.particles.X,Dp_inv);}

    example.gather_scatter.template Gather<HELPER>(true,
        [](int p,HELPER& h){h=HELPER();},
        [this,dt,&particles](int p,const PARTICLE_GRID_ITERATOR<TV>& it,HELPER& h)
        {
            T w=it.Weight();
            TV dw=it.Gradient();
            TV_INT index=it.Index();
            TV V_new=example.velocity(index);
            TV V_old=example.velocity_save(index);
            h.V_pic+=w*V_new;
            h.V_weight_old+=w*V_old;
            h.grad_Vp+=Outer_Product(V_new,dw);
            if(example.use_affine)
                h.B+=w*Outer_Product(V_new,particles.X(p)-example.grid.Center(index));
        },
        [this,dt,&particles,&Dp_inv](int p,HELPER& h)
        {
            MATRIX<T,TV::m> A=dt*h.grad_Vp+1;
            if(example.quad_F_coeff) A+=sqr(dt*h.grad_Vp)*example.quad_F_coeff;
            particles.F(p)=A*particles.F(p);

            T J=particles.F(p).Determinant();
            particles.F(p)=MATRIX<T,TV::m>()+pow<1,TV::m>(J);

            TV xp_new=particles.X(p)+dt*h.V_pic;
            if(example.use_affine){
                if(example.weights->use_gradient_transfer) particles.C(p)=h.grad_Vp;
                else particles.C(p)=h.B*Dp_inv(p);}
            particles.X(p)=xp_new;
            TV V_flip=particles.V(p)+h.V_pic-h.V_weight_old;
            particles.V(p)=V_flip*example.flip+h.V_pic*(1-example.flip);

            Reflect_Or_Invalidate_Particle(p);
        });
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
template<class TV> void MPM_MICROPOLAR_DRIVER<TV>::
Grid_To_Particle_Limit_Dt()
{
    if(!example.use_strong_cfl) return;
    struct HELPER
    {
        TV V_pic,V_pic_s,V_weight_old;
        MATRIX<T,TV::m> grad_Vp,grad_Vp_s;
    };

    T dt=example.dt,s=1;

    // TODO: this is NOT threadsafe.
    example.gather_scatter.template Gather<HELPER>(true,
        [](int p,HELPER& h){h=HELPER();},
        [this](int p,const PARTICLE_GRID_ITERATOR<TV>& it,HELPER& h)
        {
            T w=it.Weight();
            TV dw=it.Gradient();
            TV_INT index=it.Index();
            TV V_new=example.velocity(index);
            TV V_old=example.velocity_save(index);

            h.V_pic+=w*V_old;
            h.V_pic_s+=w*(V_new-V_old);
            h.V_weight_old+=w*V_old;
            TV V_grid=V_old;
            TV V_grid_s=V_new-V_old;
            h.grad_Vp+=Outer_Product(V_grid,dw);
            h.grad_Vp_s+=Outer_Product(V_grid_s,dw);
        },
        [this,dt,&s](int p,HELPER& h)
        {
            Enforce_Limit_Max(s,example.cfl_F,dt*h.grad_Vp,dt*h.grad_Vp_s);
            TV xp_new_s=dt*h.V_pic,xp_new_s2=dt*h.V_pic_s;
            Enforce_Limit_Max(s,example.cfl,xp_new_s,xp_new_s2);
        });
    if(example.dt*s<example.min_dt) s=example.min_dt/example.dt;
    if(s>=1) return;
    LOG::printf("X J CFL scale: %g -> %g\n",example.dt,example.dt*s);
    example.dt*=s;
    example.velocity.array=(example.velocity.array-example.velocity_save.array)*s+example.velocity_save.array;
}
//#####################################################################
// Function Limit_Dt_Sound_Speed
//#####################################################################
template<class TV> void MPM_MICROPOLAR_DRIVER<TV>::
Limit_Dt_Sound_Speed()
{
    if(!example.use_sound_speed_cfl) return;
    T dt=example.dt;
    LOG::printf("forces %i\n",example.forces.m);
    for(int f=0;f<example.forces.m;f++){
        if(const MPM_FINITE_ELEMENTS<TV>* force=dynamic_cast<MPM_FINITE_ELEMENTS<TV>*>(example.forces(f))){
            LOG::printf("found force\n");
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
    LOG::printf("SOUND CFL %g %g\n",example.dt,dt);
    if(dt>=example.dt) return;
    T s=dt/example.dt;
    LOG::printf("SOUND CFL %g %g (%g)\n",example.dt,dt,s);
    example.velocity.array=(example.velocity.array-example.velocity_save.array)*s+example.velocity_save.array;
    example.dt=dt;
}
//#####################################################################
// Function Apply_Forces
//#####################################################################
template<class TV> void MPM_MICROPOLAR_DRIVER<TV>::
Apply_Forces()
{
    example.Capture_Stress();
    LOG::printf("max velocity: %.16P\n",Max_Particle_Speed());
    example.Precompute_Forces(example.time,example.dt,0);
    forces.Fill(TV());
    example.Add_Forces(forces,example.time);
    Reflection_Boundary_Condition(forces,true);
    for(int i=0;i<example.valid_grid_indices.m;i++){
        int p=example.valid_grid_indices(i);
        example.velocity.array(p)+=example.dt/example.mass.array(p)*forces.array(p);}
}
//#####################################################################
// Function Compute_Dt
//#####################################################################
template<class TV> typename TV::SCALAR MPM_MICROPOLAR_DRIVER<TV>::
Compute_Dt() const
{
    T critical_speed=example.cfl*example.grid.dX.Min()/example.max_dt;
    T v=Grid_V_Upper_Bound();
    return (v>critical_speed)?(example.cfl*example.grid.dX.Min()/v):example.max_dt;
}
//#####################################################################
// Function Max_Particle_Speed
//#####################################################################
template<class TV> typename TV::SCALAR MPM_MICROPOLAR_DRIVER<TV>::
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
template<class TV> typename TV::SCALAR MPM_MICROPOLAR_DRIVER<TV>::
Grid_V_Upper_Bound() const
{
    if(!example.use_affine) return Max_Particle_Speed();
    T result=0;
    T xi=(T)6*sqrt((T)TV::m)*example.grid.one_over_dX.Min();
#pragma omp parallel for reduction(max:result)
    for(int k=0;k<example.simulated_particles.m;k++){
        int p=example.simulated_particles(k);
        T v=example.particles.V(p).Magnitude();
        if(example.particles.store_C) v+=example.particles.C(p).Frobenius_Norm()*xi;
        result=max(result,v);}
    return result;
}
//#####################################################################
// Function Update_Simulated_Particles
//#####################################################################
template<class TV> void MPM_MICROPOLAR_DRIVER<TV>::
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
}
template<class T,int d> void flip(VECTOR<T,d>& u,int a){u(a)=-u(a);}
template<class T> void flip(T& u,int a){}
//#####################################################################
// Function Reflection_Boundary_Condition
//#####################################################################
template<class TV> template<class S> void MPM_MICROPOLAR_DRIVER<TV>::
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
template<class TV> void MPM_MICROPOLAR_DRIVER<TV>::
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
namespace PhysBAM{
template class MPM_MICROPOLAR_DRIVER<VECTOR<float,2> >;
template class MPM_MICROPOLAR_DRIVER<VECTOR<float,3> >;
template class MPM_MICROPOLAR_DRIVER<VECTOR<double,2> >;
template class MPM_MICROPOLAR_DRIVER<VECTOR<double,3> >;
}
