//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Log/DEBUG_SUBSTEPS.h>
#include <Core/Log/FINE_TIMER.h>
#include <Core/Log/LOG.h>
#include <Core/Log/SCOPE.h>
#include <Core/Math_Tools/pow.h>
#include <Core/Math_Tools/RANGE_ITERATOR.h>
#include <Core/Matrices/FRAME.h>
#include <Core/Random_Numbers/RANDOM_NUMBERS.h>
#include <Tools/Krylov_Solvers/CONJUGATE_GRADIENT.h>
#include <Tools/Krylov_Solvers/MINRES.h>
#include <Tools/Nonlinear_Equations/NEWTONS_METHOD.h>
#include <Tools/Parallel_Computation/APPEND_HOLDER.h>
#include <Tools/Polynomials/QUADRATIC.h>
#include <Grid_Tools/Grids/CELL_ITERATOR.h>
#include <Grid_Tools/Grids/CELL_ITERATOR_THREADED.h>
#include <Grid_Tools/Grids/FACE_ITERATOR.h>
#include <Geometry/Seeding/POISSON_DISK_SURFACE.h>
#include <Deformables/Collisions_And_Interactions/IMPLICIT_OBJECT_COLLISION_PENALTY_FORCES.h>
#include <Deformables/Constitutive_Models/DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE.h>
#include <Deformables/Constitutive_Models/ISOTROPIC_CONSTITUTIVE_MODEL.h>
#include <Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <Deformables/Forces/FINITE_VOLUME.h>
#include <Solids/Solids/SOLID_BODY_COLLECTION.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_DRIVER.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_EXAMPLE.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_PARTICLES.h>
#include <Hybrid_Methods/Forces/MPM_FINITE_ELEMENTS.h>
#include <Hybrid_Methods/Forces/MPM_FORCE_HELPER.h>
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
template<class TV> MPM_DRIVER<TV>::
MPM_DRIVER(MPM_EXAMPLE<TV>& example)
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
template<class TV> MPM_DRIVER<TV>::
~MPM_DRIVER()
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
template<class TV> void MPM_DRIVER<TV>::
Execute_Main_Program()
{
    Step([=](){Initialize();},"initialize");
    Simulate_To_Frame(example.last_frame);
    FINE_TIMER::Dump_Timing_Info();
}
//#####################################################################
// Initialize
//#####################################################################
template<class TV> void MPM_DRIVER<TV>::
Initialize()
{
    TIMER_SCOPE_FUNC;
    LOG::cout<<std::setprecision(16)<<std::endl;
    DEBUG_SUBSTEPS::write_substeps_level=example.substeps_delay_frame<0?example.write_substeps_level:-1;

    if(example.auto_restart){
        example.viewer_dir.Read_Last_Frame(0);
        example.restart=example.viewer_dir.frame_stack(0);}
    else if(example.restart)
        example.viewer_dir.Set(example.restart);
    if(example.restart) example.viewer_dir.Make_Common_Directory(true);
    current_frame=example.restart;
    example.time=current_frame*example.frame_dt;

    example.Initialize();
    PHYSBAM_ASSERT(example.grid.Is_MAC_Grid());
    if(example.restart)
        example.Read_Output_Files();

    example.mass.Resize(example.grid.Domain_Indices(example.ghost));
    example.velocity.Resize(example.grid.Domain_Indices(example.ghost));
    example.velocity_save.Resize(example.grid.Domain_Indices(example.ghost));
    example.velocity_friction_save.Resize(example.grid.Domain_Indices(example.ghost));
    dv.u.Resize(example.grid.Domain_Indices(example.ghost));
    rhs.u.Resize(example.grid.Domain_Indices(example.ghost));
    objective.system.tmp.u.Resize(example.grid.Domain_Indices(example.ghost));

    example.particles.Store_B(example.use_affine);
    example.particles.Store_S(example.use_oldroyd);
    if(example.particles.store_B) example.Dp_inv.Resize(example.particles.X.m);

    RANGE<TV_INT> range(example.grid.Cell_Indices(example.ghost));
    example.location.Resize(range,no_init);
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

    if(!example.restart) Write_Output_Files();
    PHYSBAM_DEBUG_WRITE_SUBSTEP("after init",1);
}
//#####################################################################
// Function Advance_One_Time_Step
//#####################################################################
template<class TV> void MPM_DRIVER<TV>::
Advance_One_Time_Step()
{
    TIMER_SCOPE_FUNC;
    Step([=](){Update_Simulated_Particles();},"simulated-particles");
    Print_Particle_Stats("particle state",example.dt);
    Step([=](){Update_Particle_Weights();},"update-weights",false);
    Step([=](){example.gather_scatter.Prepare_Scatter(example.particles);},"prepare-scatter",false);
    Step([=](){Particle_To_Grid();},"p2g");
    Print_Grid_Stats("after particle to grid",example.dt,example.velocity,0);
    Print_Energy_Stats("after particle to grid",example.velocity);
    Step([=](){Apply_Forces();},"forces");
    Print_Grid_Stats("after forces",example.dt,example.velocity,&example.velocity_save);
    Step([=](){Reduce_Dt();},"reduce-dt",false);
    Step([=](){Grid_To_Particle();},"g2p",true);
    Step([=](){Update_Plasticity_And_Hardening();},"plasticity-hardening",false);
}
//#####################################################################
// Function Simulate_To_Frame
//#####################################################################
template<class TV> void MPM_DRIVER<TV>::
Simulate_To_Frame(const int frame)
{
    TIMER_SCOPE_FUNC;
    for(;current_frame<frame;current_frame++){
        LOG::SCOPE scope("FRAME","frame %d",current_frame+1);
        for(int i=0;i<example.begin_frame.m;i++)
            example.begin_frame(i)(current_frame);
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
            T original_dt=example.dt;
            LOG::cout<<"substep dt: "<<example.dt<<std::endl;
            PHYSBAM_ASSERT(example.dt>0);

            Step([=](){Advance_One_Time_Step();},"time-step");
            LOG::cout<<"actual dt: "<<example.dt<<std::endl;

            // Time step was reduced
            if(example.dt<original_dt){
                LOG::printf("dt reduced: %g to %g  (%g)\n",original_dt,example.dt,original_dt-example.dt);
                next_time=example.time+example.dt;
                done=false;}

            PHYSBAM_DEBUG_WRITE_SUBSTEP("end substep %i",0,substep);
            example.time=next_time;}
        for(int i=0;i<example.end_frame.m;i++)
            example.end_frame(i)(current_frame);
        Write_Output_Files();}
}
//#####################################################################
// Function Write_Substep
//#####################################################################
template<class TV> void MPM_DRIVER<TV>::
Write_Substep(const std::string& title)
{
    TIMER_SCOPE_FUNC;
    example.frame_title=title;
    LOG::printf("Writing substep [%s]: output_number=%P, time=%g, frame=%i\n",
        example.frame_title,example.viewer_dir.frame_stack,example.time,current_frame);
    Write_Output_Files();
    example.frame_title="";
}
//#####################################################################
// Write_Output_Files
//#####################################################################
template<class TV> void MPM_DRIVER<TV>::
Write_Output_Files()
{
    TIMER_SCOPE_FUNC;
    LOG::SCOPE scope("Write_Output_Files");
    example.viewer_dir.Start_Directory(0,example.frame_title);
    example.frame_title="";
    example.Write_Output_Files();
    if(example.extra_render && example.r_F){
        if (ARRAY_VIEW<T>* prop4r=example.particles.template Get_Array<T>("prop4r")) (*prop4r).Fill((T)1);}
    example.viewer_dir.Finish_Directory();
}
//#####################################################################
// Function Update_Particle_Weights
//#####################################################################
template<class TV> void MPM_DRIVER<TV>::
Update_Particle_Weights()
{
    TIMER_SCOPE_FUNC;
    example.weights->Update(example.particles.X);
    if(example.particles.store_B){
        example.Dp_inv.Resize(example.particles.X.m);
        example.weights->Dp_Inverse(example.particles.X,example.Dp_inv);}
}
//#####################################################################
// Function Particle_To_Grid
//#####################################################################
template<class TV> void MPM_DRIVER<TV>::
Particle_To_Grid()
{
    TIMER_SCOPE_FUNC;
    MPM_PARTICLES<TV>& particles=example.particles;

#pragma omp parallel for
    for(int i=0;i<example.mass.array.m;i++){
        example.mass.array(i)=0;
        example.velocity.array(i)=TV();
        example.velocity_save.array(i)=TV();
        example.velocity_friction_save.array(i)=TV();}

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
    PHYSBAM_DEBUG_WRITE_SUBSTEP("momentum",1);
    if(example.reflection_bc && !example.use_full_reflection)
    {
        Reflect_Boundary_Mass_Momentum();
        PHYSBAM_DEBUG_WRITE_SUBSTEP("after reflect m mv",1);
    }

    APPEND_HOLDER<int> flat_h(example.valid_grid_indices);
    APPEND_HOLDER<TV_INT> indices_h(example.valid_grid_cell_indices);
#pragma omp parallel
    {
#pragma omp single
        {
            flat_h.Init();
            indices_h.Init();
        }
#pragma omp barrier
        ARRAY<int>& flat_t=flat_h.Array();
        ARRAY<TV_INT>& indices_t=indices_h.Array();
        for(CELL_ITERATOR_THREADED<TV> it(example.grid,example.ghost);it.Valid();it.Next()){
            int i=example.mass.Standard_Index(it.index);
            if(example.mass.array(i)){
                flat_t.Append(i);
                indices_t.Append(it.index);
                TV v=example.velocity.array(i)/example.mass.array(i);
                example.velocity.array(i)=v;
                example.velocity_save.array(i)=v;}
            else example.velocity.array(i)=TV();}
    }
    flat_h.Combine();
    indices_h.Combine();
}
//#####################################################################
// Function Reflection_Boundary_Mass_Momentum
//#####################################################################
template <class TV> void MPM_DRIVER<TV>::
Reflect_Boundary_Mass_Momentum()
{
    TIMER_SCOPE_FUNC;
    TV_INT ranges[2]={TV_INT(),example.grid.numbers_of_cells};
    for(RANGE_ITERATOR<TV::m> it(example.grid.Domain_Indices(),example.ghost,0,
            RI::ghost|RI::duplicate_corners);it.Valid();it.Next()){
        auto bc_type=example.side_bc_type(it.side);
        TV_INT image=it.index;
        int axis=it.side/2,side=it.side%2;
        image(axis)=2*ranges[side](axis)-it.index(axis)-1;
        T& m_i=example.mass(image),&m_o=example.mass(it.index),m_t=m_i+m_o;
        TV& mv_i=example.velocity(image),&mv_o=example.velocity(it.index);
        TV mv_id=mv_i-mv_o,mv_od=-mv_id,mv_in=mv_i+mv_o,mv_on=mv_in;
        if(bc_type!=example.BC_FREE && example.bc_velocity){
            TV nearest_point=example.grid.Clamp(example.grid.Center(it.index));
            TV bc=example.bc_velocity(nearest_point,example.time);
            mv_id+=2*m_o*bc;
            mv_od+=2*m_i*bc;}
        if(bc_type==example.BC_SLIP){
            mv_in(axis)=mv_id(axis);
            mv_on(axis)=mv_od(axis);}
        if(bc_type==example.BC_NOSLIP){
            mv_i=mv_id;
            mv_o=mv_od;}
        else{
            mv_i=mv_in;
            mv_o=mv_on;}
        m_o=m_i=m_t;};
}
template<class TV,class T>
PAIR<TV,TV> Reflect_Helper(T m_i,T m_o,TV v_i,TV v_o,auto type,T mu,int a,int s,TV bc)
{
    T m_t=m_i+m_o,a_i=m_i/m_t,a_o=1-a_i;
    TV v_in=a_i*v_i+a_o*v_o,v_on=v_in;
    if(type==MPM_EXAMPLE<TV>::BC_FREE) return {v_in,v_on};

    TV v_id=a_i*v_i+a_o*(2*bc-v_o),v_od=a_o*v_o+a_i*(2*bc-v_i);
    if(type==MPM_EXAMPLE<TV>::BC_NOSLIP) return {v_id,v_od};

    if(type==MPM_EXAMPLE<TV>::BC_SLIP){
        v_in(a)=v_id(a);
        v_on(a)=v_od(a);
        return {v_in,v_on};}

    T d_v=2*a_o*a_i*(2*bc(a)-v_i(a)-v_o(a));
    T d_v_n=s?-d_v:d_v;
    if(d_v_n<0) return {v_i,v_o};

    TV v_id_n,v_od_n;
    v_id_n(a)=v_id(a);
    v_od_n(a)=v_od(a);
    v_in(a)=0;

    T d_v_t_mag=v_in.Normalize();
    T r=d_v_t_mag-mu*d_v_n;
    if(r<0) return {v_id_n,v_od_n};
    return {v_id_n+r*v_in,v_od_n+r*v_in};
}

//#####################################################################
// Function Reflect_Boundary_Friction
//#####################################################################
template <class TV> void MPM_DRIVER<TV>::
Reflect_Boundary_Friction(const ARRAY<T,TV_INT>& mass,ARRAY<TV,TV_INT>& u) const
{
    TIMER_SCOPE_FUNC;
    TV_INT ranges[2]={TV_INT(),example.grid.numbers_of_cells};
    for(RANGE_ITERATOR<TV::m> it(example.grid.Domain_Indices(),example.ghost,0,
            RI::ghost|RI::duplicate_corners);it.Valid();it.Next()){
        auto bc_type=example.side_bc_type(it.side);
        TV_INT image=it.index;
        int axis=it.side/2,side=it.side%2;
        image(axis)=2*ranges[side](axis)-it.index(axis)-1;
        T m_i=example.mass(image),m_o=example.mass(it.index);
        if(!m_i || !m_o) continue;
        TV& v_i=example.velocity(image),&v_o=example.velocity(it.index);
        TV bc;
        if(bc_type!=example.BC_FREE && example.bc_velocity)
            bc=example.bc_velocity(example.grid.Clamp(example.grid.Center(it.index)),example.time);
        auto pr=Reflect_Helper(m_i,m_o,v_i,v_o,bc_type,example.reflection_bc_friction,axis,side,bc);
        v_i=pr.x;
        v_o=pr.y;}
}
//#####################################################################
// Function Reflect_Boundary_Velocity
//#####################################################################
template <class TV> void MPM_DRIVER<TV>::
Reflect_Boundary_Velocity(ARRAY<TV,TV_INT>& u)
{
    TIMER_SCOPE_FUNC;
    if(!example.reflection_bc) return;
    TV_INT ranges[2]={TV_INT(),example.grid.numbers_of_cells};
    for(RANGE_ITERATOR<TV::m> it(example.grid.Domain_Indices(),example.ghost,0,
            RI::ghost|RI::duplicate_corners);it.Valid();it.Next()){
        auto bc_type=example.side_bc_type(it.side);
        TV_INT image=it.index;
        int axis=it.side/2,side=it.side%2;
        image(axis)=2*ranges[side](axis)-it.index(axis)-1;
        T m=example.mass(it.index);
        if(!m) continue;
        TV v_i=u(image),&v_o=u(it.index);
        TV v_od=-v_i,v_on=v_i;
        if(bc_type!=example.BC_FREE && example.bc_velocity){
            TV nearest_point=example.grid.Clamp(example.grid.Center(it.index));
            TV bc=example.bc_velocity(nearest_point,example.time);
            v_od+=bc*2;}
        if(bc_type==example.BC_SLIP) v_on(axis)=v_od(axis);
        if(bc_type==example.BC_NOSLIP) v_o=v_od;
        else v_o=v_on;}
}
//#####################################################################
// Function Reflect_Boundary_Force
//#####################################################################
template <class TV> void MPM_DRIVER<TV>::
Reflect_Boundary_Force(ARRAY<TV,TV_INT>& force)
{
    TIMER_SCOPE_FUNC;
    if(!example.reflection_bc) return;
    TV_INT ranges[2]={TV_INT(),example.grid.numbers_of_cells};
    for(RANGE_ITERATOR<TV::m> it(example.grid.Domain_Indices(),example.ghost,0,RI::ghost|RI::duplicate_corners);it.Valid();it.Next()){
        auto bc_type=example.side_bc_type(it.side);
        TV_INT image=it.index;
        int axis=it.side/2,side=it.side%2;
        image(axis)=2*ranges[side](axis)-it.index(axis)-1;
        TV& f_i=force(image),&f_o=force(it.index);
        TV f_id=f_i-f_o,f_od=-f_id,f_in=f_i+f_o,f_on=f_in;
        if(bc_type==example.BC_SLIP){
            f_in(axis)=f_id(axis);
            f_on(axis)=f_od(axis);}
        if(bc_type==example.BC_NOSLIP){
            f_i=f_id;
            f_o=f_od;}
        else{
            f_i=f_in;
            f_o=f_on;}}
}
//#####################################################################
// Function Grid_To_Particle
//#####################################################################
template<class TV> void MPM_DRIVER<TV>::
Grid_To_Particle()
{
    TIMER_SCOPE_FUNC;
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
            if(example.dilation_only)
            {
                T J=particles.F(p).Determinant()*(dt*h.grad_Vp.Trace()+1);
                particles.F(p)=MATRIX<T,TV::m>()+pow<1,TV::m>(J);
            }
            else particles.F(p)=A*particles.F(p);
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
}
//#####################################################################
// Function Reduce_Dt
//#####################################################################
template<class TV> void MPM_DRIVER<TV>::
Reduce_Dt()
{
    TIMER_SCOPE_FUNC;
    T dt=example.dt;
    if(example.use_strong_cfl)
    {
        T new_dt=Grid_To_Particle_Limit_Dt();
        if(new_dt<dt)
        {
            LOG::printf("STRONG CFL %g %g (%g)\n",new_dt,dt,new_dt/dt);
            dt=new_dt;
        }
    }
    dt=std::min(dt,Limit_Dt_Sound_Speed());
    dt=std::max(dt,example.min_dt);
    if(dt>=example.dt*(T)0.99) return;

    T s=dt/example.dt;
    LOG::printf("dt scale: %g -> %g (%g)\n",example.dt,dt,s);
#pragma omp parallel for
    for(int i=0;i<example.valid_grid_indices.m;i++){
        int j=example.valid_grid_indices(i);
        TV vs=example.velocity_save.array(j);
        TV& v=example.velocity.array(j);
        v=(v-vs)*s+vs;
        TV& vf=example.velocity_friction_save.array(j);
        vf=(vf-vs)*s+vs;}
    example.dt=dt;
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
template<class TV> auto MPM_DRIVER<TV>::
Grid_To_Particle_Limit_Dt() -> T
{
    TIMER_SCOPE_FUNC;
    struct HELPER
    {
        TV V_pic,V_pic_s,V_weight_old;
        MATRIX<T,TV::m> grad_Vp,grad_Vp_s;
        T s=100;
    };

    T dt=example.dt;
    T midpoint_frac=example.use_midpoint?(T).5:1;
    T s=100;
    
    example.gather_scatter.template Gather<HELPER>(true,
        [&s](HELPER& h){s=std::min(s,h.s);},
        [](int p,HELPER& h){T x=h.s;h=HELPER();h.s=x;},
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
        [this,dt](int p,HELPER& h)
        {
            T s_save=h.s;
            h.s=100;
            if(example.dilation_only) Enforce_Limit_Max(h.s,example.cfl_F,dt*h.grad_Vp.Trace(),dt*h.grad_Vp_s.Trace());
            else Enforce_Limit_Max(h.s,example.cfl_F,dt*h.grad_Vp,dt*h.grad_Vp_s);
            TV xp_new_s,xp_new_s2;
            if(example.use_midpoint){xp_new_s=dt/2*(h.V_weight_old+h.V_pic);xp_new_s2=dt/2*h.V_pic_s;}
            else{xp_new_s=dt*h.V_pic;xp_new_s2=dt*h.V_pic_s;}
            Enforce_Limit_Max(h.s,example.cfl,xp_new_s,xp_new_s2);
            if(example.extra_render && example.r_F){
                if (ARRAY_VIEW<T>* prop4r=example.particles.template Get_Array<T>("prop4r")) {
                    if((T)0==(*prop4r)(p)){
                        (*prop4r)(p)=h.s;}
                    else if((*prop4r)(p)>h.s) (*prop4r)(p)=h.s;}}
            h.s=min(h.s,s_save);
        });
    if(example.verbose_cfl) LOG::printf("F CFL %P\n",example.dt*s);
    return example.dt*s;
}
//#####################################################################
// Function Limit_Dt_Sound_Speed
//#####################################################################
template<class TV> auto MPM_DRIVER<TV>::
Limit_Dt_Sound_Speed() -> T
{
    TIMER_SCOPE_FUNC;
    T dt=example.dt;
    if(example.use_sound_speed_cfl)
    {
        T max_speed=Compute_Max_Sound_Speed();
        LOG::printf("max sound speed: %.16P\n",max_speed);
        LOG::printf("dx: %.16P\n",example.grid.dX.Min());
        LOG::printf("dx/soundspeed: %.16P\n",Robust_Divide(example.grid.dX.Min(),max_speed));
        T new_dt=Robust_Divide(example.grid.dX.Min(),max_speed)*example.cfl_sound;
        LOG::printf("SOUND CFL %g %g (%g)\n",example.dt,new_dt,Robust_Divide(new_dt,example.dt));
        dt=std::min(dt,new_dt);
    }
    if(example.use_single_particle_cfl && example.dilation_only)
    {
        T new_dt=Max_Dt_Single_Particle_Pressure();
        if(new_dt<dt)
        {
            dt=new_dt;
            LOG::printf("SINGLE PARTICLE (J) %g %g (%g)\n",example.dt,dt,dt/example.dt);
        }
    }
    if(example.use_single_particle_cfl && !example.dilation_only)
    {
        T new_dt=Max_Dt_Single_Particle();
        if(new_dt<dt)
        {
            dt=new_dt;
            LOG::printf("SINGLE PARTICLE (F) %g %g (%g)\n",example.dt,dt,dt/example.dt);
        }
    }
    return dt;
}
//#####################################################################
// Function Max_Dt_Single_Particle_Pressure
//#####################################################################
template<class TV> auto MPM_DRIVER<TV>::
Max_Dt_Single_Particle_Pressure() const -> T
{
    TIMER_SCOPE_FUNC;
    T dt=example.dt;
    T K=TV::m==2?6:(T)3.14;
    for(int f=0;f<example.forces.m;f++){
        if(const MPM_FINITE_ELEMENTS<TV>* force=dynamic_cast<MPM_FINITE_ELEMENTS<TV>*>(example.forces(f))){
#pragma omp parallel for reduction(min:dt)
            for(int k=0;k<example.simulated_particles.m;k++){
                int p=example.simulated_particles(k);
                T density=example.particles.mass(p)/example.particles.volume(p);
                T J=force->sigma(p).Determinant();
                T num=sqr(example.cfl_single_particle*example.grid.dX.Min())*density;
                T den=K*TV::m;
#if 0
                T dp_over_Jm1=force->constitutive_model.Robust_Divided_Pressure(J,p);
                den*=sqr(J)*dp_over_Jm1*J/(J+1);
#else
                T lambda=force->constitutive_model.Pressure_Bound(J,p);
                if(J>=1)
                {
                    num*=J+1;
                    den*=J*J*J*lambda;
                }
                else
                {
                    num*=2;
                    den*=sqr(2-J)*lambda;
                }
#endif
                if(num<sqr(dt)*den) dt=sqrt(num/den);}}}
    PHYSBAM_ASSERT(example.lagrangian_forces.m==0);
    return dt;
}
//#####################################################################
// Function Max_Dt_Single_Particle
//#####################################################################
template<class TV> auto MPM_DRIVER<TV>::
Max_Dt_Single_Particle() const -> T
{
    PHYSBAM_FATAL_ERROR();
    return 0;
}
//#####################################################################
// Function Conjugate_Stress_Diff
//#####################################################################
template<class TV> SYMMETRIC_MATRIX<typename TV::SCALAR,TV::m>
Conjugate_Stress_Diff(const DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<TV>& dpdf,const TV& u)
{
    typedef typename TV::SCALAR T;
    DIAGONAL_MATRIX<T,TV::m> U(u);
    SYMMETRIC_MATRIX<T,TV::m> S=SYMMETRIC_MATRIX<T,TV::m>::Conjugate(U,dpdf.H);
    for(int i=0;i<TV::SPIN::m;i++){
        int j=(i+1)%TV::m,k=(j+1)%TV::m;
        S(j,j)+=dpdf.B(i)*u(k)*u(k);
        S(k,k)+=dpdf.B(i)*u(j)*u(j);
        S(j,k)+=dpdf.C(i)*u(j)*u(k);}
    return S;
}
template<class TV,class T>
inline T Compute_Maximum_Tensor_Contraction_Brute(const DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<TV>& dpdf,const DIAGONAL_MATRIX<T,TV::m>& sigma)
{
    T epsilon=std::numeric_limits<T>::epsilon();
    T tolerance=sqrt(epsilon)/100;
    RANDOM_NUMBERS<T> random;
    TV u,v;
    T eig_est=-FLT_MAX;
    for(int i=0;i<100;i++){
        TV u0,v0;
        random.Fill_Uniform(u0,-1,1);
        random.Fill_Uniform(v0,-1,1);
        u0.Normalize();
        v0.Normalize();
        T e=u0.Dot(Conjugate_Stress_Diff(dpdf,sigma*v0)*u0);
        if(e>eig_est){
            eig_est=e;
            u=u0;
            v=v0;}}

    T ev=FLT_MAX,pev=0;
    int iter=10000;
    while(abs(ev-pev)>tolerance*ev && abs(ev)>tolerance){
        if(!--iter)
        {
            break;
            LOG::printf("STOP EARLY\n");
            return ev;
        }
        DIAGONAL_MATRIX<T,TV::m> E;
        MATRIX<T,TV::m> V;
        auto A=Conjugate_Stress_Diff(dpdf,u);
        A=A.Conjugate(sigma,A);
        A.Fast_Solve_Eigenproblem(E,V);
        v=V.Column(E.x.Arg_Max());
        auto B=Conjugate_Stress_Diff(dpdf,sigma*v);
        B.Fast_Solve_Eigenproblem(E,V);
        const int kk=E.x.Arg_Max();
        u=V.Column(kk);
        pev=ev;
        ev=E.x(kk);}

    T best_a=-FLT_MAX;
    for(int i=0;i<TV::m;i++)
        for(int j=0;j<TV::m;j++)
        {
            TV u0,v0;
            u0(i)=1;
            v0(j)=1;
            T e=u0.Dot(Conjugate_Stress_Diff(dpdf,sigma*v0)*u0);
            if(e>best_a)
                best_a=e;
        }
    if(best_a*1.0000001>=ev) return best_a;

#pragma omp critical
    {
        LOG::printf("XX %P %P -> %P   %P\n",u,v,ev,Robust_Divide(ev,best_a));
    }
    return ev;
}
template<class T>
inline T Compute_Maximum_Tensor_Contraction(const SYMMETRIC_MATRIX<T,2>& H, const VECTOR<T,1>& A0, const VECTOR<T,1>& A1, const VECTOR<T,1>& B)
{
    T A01=A0.x,A10=A1.x,H00=H.x00,H10=H.x10+B.x,H11=H.x11;
    T t0=A01*A10-A01*H00-A10*H11+H00*H11-H10*H10;
    T t1=-H11+A01-2*H10-H00+A10;
    T t2=-H11+A01+2*H10-H00+A10;
    T t4=A10*A10-A10*H00-A10*H11+H00*H11-H10*H10;
    T t5=A01*A10-A01*H11-A10*H11-H10*H10+H11*H11;
    T t6=A01*A10-A01*H00-A10*H00+H00*H00-H10*H10;
    QUADRATIC<T> q(t0*t1*t2,-2*(t4+t5)*t0,t4*t5);
    q.Compute_Roots();

    T x = max(H.x00,H.x11,A0.x,A1.x);
    return x;
    auto do_v02=[=,&x](T v02)
    {
        T v12=1-v02;
        T den=(t5*v12-t6*v02);
        T num=v02*(2*t0*v12-t6);
        if(den<0)
        {
            den=-den;
            num=-num;
        }

        if(num>=0 && num<=den)
        {
            T nu=(num-v12*den)*H10;
            T de=((A01-H00)*num+(H00-A10)*v12*den);
            if(abs(de)>1e-10){
                T t0=nu/de;
                T e=(((A10*v12/v02+H00)*t0+2*H10)*v12*t0+A01*v02+H11*v12)*num/den;
                if(e>x)
                {
                    x=e;
                    LOG::printf("SPECIAL\n");
                }}
            else{
                LOG::printf("Not Stable!\n");}
        }
    };
    T epsilon=std::numeric_limits<T>::epsilon()/100;
    for(int i=0;i<2;i++)
        if(q.roots>i && q.root[i]>=epsilon && q.root[i]<=1-epsilon)
            do_v02(q.root[i]);
    return x;
}
template<class TV,class T>
inline T Compute_Maximum_Tensor_Contraction(const DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<TV>& dpdf,const DIAGONAL_MATRIX<T,2>& sigma)
{
    auto S=SYMMETRIC_MATRIX<T,TV::m>::Conjugate(sigma,dpdf.H);
    VECTOR<T,1> B0(dpdf.B(0)*sigma(1)*sigma(1)),B1(dpdf.B(0)*sigma(0)*sigma(0)),C(dpdf.C(0)*sigma(1)*sigma(0));
    return Compute_Maximum_Tensor_Contraction(S,B0,B1,C);
}
template<class T,class TV>
const T Compute_Maximum_Tensor_Contraction(const DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<TV>& dpdf,const DIAGONAL_MATRIX<T,3>& sigma)
{
    DIAGONAL_MATRIX<T,2> S2;
    SYMMETRIC_MATRIX<T,2> H2;
    T ev=-FLT_MAX;
    for(int d=0;d<TV::m;++d){
        for(int i=0;i<2;++i){
            S2(i)=sigma((d+i+1)%TV::m);}
        H2.x00=dpdf.H((d+1)%TV::m,(d+1)%TV::m);
        H2.x11=dpdf.H((d+2)%TV::m,(d+2)%TV::m);
        H2.x10=dpdf.H((d+1)%TV::m,(d+2)%TV::m);
        auto S=SYMMETRIC_MATRIX<T,2>::Conjugate(S2,H2);
        VECTOR<T,1> B0(dpdf.B(d)*S2(1)*S2(1)),B1(dpdf.B(d)*S2(0)*S2(0)),C(dpdf.C(d)*S2(1)*S2(0));
        ev=max(ev,Compute_Maximum_Tensor_Contraction(S,B0,B1,C));}
    return ev;
}
//#####################################################################
// Function Test_Sound_Speed
//#####################################################################
template<class TV> void MPM_DRIVER<TV>::
Test_Sound_Speed(int num) const
{
    TIMER_SCOPE_FUNC;
    DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<TV> dpdf;
    DIAGONAL_MATRIX<T,TV::m> sigma;
    RANDOM_NUMBERS<T> r;
    LOG::printf("===== TEST %i =====\n",num);
    for(int i=0;i<num;i++)
    {
        r.Fill_Uniform(dpdf.H,-1,1);
        r.Fill_Uniform(dpdf.B,-1,1);
        r.Fill_Uniform(dpdf.C,-1,1);
        r.Fill_Uniform(sigma.x,-1,1);
        T x=Compute_Maximum_Tensor_Contraction(dpdf,sigma);
        T y=Compute_Maximum_Tensor_Contraction_Brute(dpdf,sigma);
        T r=abs(x-y)/maxabs(maxabs(x,y),(T)FLT_MIN);
        if(r>1e-6) LOG::printf("COMP: %P %P -> %P\n",x,y,r);
    }
}
//#####################################################################
// Function Compute_Max_Sound_Speed
//#####################################################################
template<class TV> auto MPM_DRIVER<TV>::
Compute_Max_Sound_Speed() const -> T
{
    TIMER_SCOPE_FUNC;
    T max_speed=0;
    for(int f=0;f<example.forces.m;f++){
        if(const MPM_FINITE_ELEMENTS<TV>* force=dynamic_cast<MPM_FINITE_ELEMENTS<TV>*>(example.forces(f))){
#pragma omp parallel for reduction(max:max_speed)
            for(int k=0;k<example.simulated_particles.m;k++){
                int p=example.simulated_particles(k);
                T density=example.particles.mass(p)/example.particles.volume(p);
                T speed=ISOTROPIC_CONSTITUTIVE_MODEL<T,TV::m>::Sound_Speed(
                    force->sigma(p),force->dPi_dF(p),density);
                    if(example.extra_render){
                        if(example.r_sound_speed){
                            if (ARRAY_VIEW<T>* prop4r=example.particles.template Get_Array<T>("prop4r")) (*prop4r)(p)=speed;}
                        if(example.r_cfl){
                            if (ARRAY_VIEW<T>* prop4r=example.particles.template Get_Array<T>("prop4r")){
                                T s=Robust_Divide(example.grid.dX.Min(),speed)/example.dt;
                                (*prop4r)(p)=s<1?s:(*prop4r)(p);}}} 
                max_speed=max(max_speed,speed);}}}
    for(auto* pf:example.lagrangian_forces)
        if(FINITE_VOLUME<TV,TV::m>* fv=dynamic_cast<FINITE_VOLUME<TV,TV::m>*>(pf))
            max_speed=max(max_speed,fv->Compute_Sound_Speed());
    return max_speed;
}
//#####################################################################
// Function Update_Plasticity_And_Hardening
//#####################################################################
template<class TV> void MPM_DRIVER<TV>::
Update_Plasticity_And_Hardening()
{
    TIMER_SCOPE_FUNC;
    for(int i=0;i<example.plasticity_models.m;i++)
        example.plasticity_models(i)->Update_Particles();
}
//#####################################################################
// Function Apply_Particle_Forces
//#####################################################################
template<class TV> void MPM_DRIVER<TV>::
Apply_Particle_Forces(ARRAY<TV,TV_INT>& F)
{
    TIMER_SCOPE_FUNC;
    if(example.lagrangian_forces.m) {
        example.lagrangian_forces_F.Resize(example.particles.X.m,no_init);
#pragma omp parallel for
        for(int i=0;i<example.lagrangian_forces_F.m;i++)
            example.lagrangian_forces_F(i)=TV();
        example.solid_body_collection.deformable_body_collection.Add_Velocity_Independent_Forces(example.lagrangian_forces_F,example.time);
        example.gather_scatter.template Scatter<int>(false,
            [this,&F](int p,const PARTICLE_GRID_ITERATOR<TV>& it,int tid){
                F(it.Index())+=it.Weight()*example.lagrangian_forces_F(p);});}
}
//#####################################################################
// Function Apply_Grid_Forces
//#####################################################################
template<class TV> void MPM_DRIVER<TV>::
Apply_Grid_Forces(ARRAY<TV,TV_INT>& F)
{
    TIMER_SCOPE_FUNC;
    for(int i=0;i<example.forces.m;i++)
        example.forces(i)->Add_Forces(F,example.time);
}
//#####################################################################
// Function Apply_Forces
//#####################################################################
template<class TV> void MPM_DRIVER<TV>::
Apply_Forces()
{
    TIMER_SCOPE_FUNC;
    example.Capture_Stress();
    objective.Reset();
    LOG::printf("max velocity: %.16P\n",Max_Particle_Speed());
    if(example.use_symplectic_euler){
        example.force_helper.B.Resize(example.particles.number);
        example.force_helper.B.Fill(MATRIX<T,TV::m>());
        example.Precompute_Forces(example.time,example.dt,0);
#pragma omp parallel for
        for(int i=0;i<objective.tmp2.u.array.m;i++)
            objective.tmp2.u.array(i)=TV();

        TIMER_SCOPE("Apply_Forces B");
        //Add Particle Forces
        auto &F=objective.tmp2.u;
        Apply_Particle_Forces(F);
        Apply_Grid_Forces(F);
        example.velocity.Exchange(F);
        PHYSBAM_DEBUG_WRITE_SUBSTEP("raw forces",1);
        example.velocity.Exchange(F);
        if(example.reflection_bc && !example.use_full_reflection){
            Reflect_Boundary_Force(F);
            example.velocity.Exchange(F);
            PHYSBAM_DEBUG_WRITE_SUBSTEP("forces after reflect",1);
            example.velocity.Exchange(F);}

        TIMER_SCOPE("Apply_Forces C");
        for(int i=0;i<example.valid_grid_indices.m;i++){
            int p=example.valid_grid_indices(i);
            dv.u.array(p)=example.dt/example.mass.array(p)*objective.tmp2.u.array(p);}
        objective.system.forced_collisions.Remove_All();
        objective.tmp1=dv;

        if(!example.use_reflection_collision_objects)
            objective.Adjust_For_Collision(dv);

        objective.tmp0=dv;
        if((example.reflection_bc && example.use_full_reflection) ||
            example.use_reflection_collision_objects)
        {
#pragma omp parallel for
            for(int i=0;i<example.valid_grid_indices.m;i++){
                int j=example.valid_grid_indices(i);
                example.velocity.array(j)=dv.u.array(j)+objective.v0.u.array(j);}
            PHYSBAM_DEBUG_WRITE_SUBSTEP("velocities before local collisions",1);
            if(example.use_reflection_collision_objects)
            {
                Apply_Reflection_Collision_Objects();
                PHYSBAM_DEBUG_WRITE_SUBSTEP("velocities after object collisions",1);
            }
            if(example.reflection_bc && example.use_full_reflection)
            {
                Reflect_Boundary_Friction(example.mass,example.velocity);
                PHYSBAM_DEBUG_WRITE_SUBSTEP("velocities after boundary collisions",1);
            }
#pragma omp parallel for
            for(int i=0;i<example.valid_grid_indices.m;i++){
                int j=example.valid_grid_indices(i);
                dv.u.array(j)=example.velocity.array(j)-objective.v0.u.array(j);}
        }
        Apply_Friction();
    }
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
        newtons_method.require_one_iteration=!objective.Initial_Guess(dv,newtons_method.tolerance,example.asymmetric_system && !newtons_method.use_gradient_magnitude_objective,newtons_method.use_gradient_magnitude_objective);
        if(example.test_diff) objective.Test_Diff(dv);

        objective.system.forced_collisions.Remove_All();

        bool converged=newtons_method.Newtons_Method(objective,objective.system,dv,av);

        if(!converged) LOG::cout<<"WARNING: Newton's method did not converge"<<std::endl;
        Apply_Friction();
        objective.Restore_F();}
    TIMER_SCOPE("Apply_Forces E");
#pragma omp parallel for
    for(int i=0;i<example.valid_grid_indices.m;i++){
        int j=example.valid_grid_indices(i);
        example.velocity.array(j)=dv.u.array(j)+objective.v0.u.array(j);
        example.velocity_friction_save.array(j)+=objective.v0.u.array(j);}
    TIMER_SCOPE("Apply_Forces F");
    example.velocity_friction_save.array.Subset(objective.system.stuck_nodes)=objective.system.stuck_velocity;
    if(example.reflection_bc && !example.use_full_reflection){
        Reflect_Boundary_Velocity(example.velocity);
        Reflect_Boundary_Velocity(example.velocity_friction_save);}
}
//#####################################################################
// Function Apply_Friction
//#####################################################################
template<class TV> void MPM_DRIVER<TV>::
Apply_Friction()
{
    TIMER_SCOPE_FUNC;
#pragma omp parallel for
    for(int i=0;i<example.valid_grid_indices.m;i++){
        int j=example.valid_grid_indices(i);
        example.velocity_friction_save.array(j)=dv.u.array(j);}
    if(!example.collision_objects.m || example.use_reflection_collision_objects) return;
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
template<class TV> typename TV::SCALAR MPM_DRIVER<TV>::
Compute_Dt() const
{
    T critical_speed=example.cfl*example.grid.dX.Min()/example.max_dt;
    T v=Grid_V_Upper_Bound();
    return (v>critical_speed)?(example.cfl*example.grid.dX.Min()/v):example.max_dt;
}
//#####################################################################
// Function Max_Particle_Speed
//#####################################################################
template<class TV> typename TV::SCALAR MPM_DRIVER<TV>::
Max_Particle_Speed() const
{
    T v2=0;
    T max_F=0;
#pragma omp parallel for reduction(max:v2)
    for(int k=0;k<example.simulated_particles.m;k++){
        int p=example.simulated_particles(k);
        v2=max(v2,example.particles.V(p).Magnitude_Squared());
        auto F=example.particles.F(p);
        max_F=max(max_F,(F.Transpose_Times(F)-1).Frobenius_Norm_Squared());}
    LOG::printf("maximum F magnitude: %P\n",sqrt(max_F));
    return sqrt(v2);
}
//#####################################################################
// Function Grid_V_Upper_Bound
//#####################################################################
template<class TV> typename TV::SCALAR MPM_DRIVER<TV>::
Grid_V_Upper_Bound() const
{
    if(!example.use_affine || !example.use_affine_cfl) return Max_Particle_Speed();
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
template<class TV> void MPM_DRIVER<TV>::
Update_Simulated_Particles()
{
    TIMER_SCOPE_FUNC;
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
}
//#####################################################################
// Function Print_Grid_Stats
//#####################################################################
template<class TV> void MPM_DRIVER<TV>::
Print_Grid_Stats(const char* str,T dt,const ARRAY<TV,TV_INT>& u,const ARRAY<TV,TV_INT>* u0)
{
    TIMER_SCOPE_FUNC;
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
Print_Particle_Stats(const char* str,T dt)
{
    TIMER_SCOPE_FUNC;
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
    TIMER_SCOPE_FUNC;
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
// Function Reflect_Particles
//#####################################################################
template<class TV> void MPM_DRIVER<TV>::
Reflect_Or_Invalidate_Particle(int p)
{
    int f=0;
    for(int i=0;i<2*TV::m;i++)
        if(example.side_bc_type(i)!=example.BC_FREE)
            f|=1<<i;
    TV A=example.grid.domain.min_corner;
    TV B=example.grid.domain.max_corner;
    TV eps=std::numeric_limits<T>::epsilon()*4*(B-A);
    A+=eps;
    B-=eps; // avoid being on the edge and rasterizing to a grid cell outside the domain.
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
// Function Apply_Reflection_Collision_Objects
//#####################################################################
template<class TV> void MPM_DRIVER<TV>::
Apply_Reflection_Collision_Objects()
{
    if(!example.use_reflection_collision_objects) return;
    if(!example.colliding_nodes.array.m)
        example.colliding_nodes.Resize(example.mass.domain,use_init,-1);

    for(int i=example.collision_objects.m;i<example.collision_objects_reflection.m;i++)
        delete example.collision_objects_reflection(i);
    example.collision_objects_reflection.Resize(example.collision_objects.m);
    struct NODE_DATA
    {
        int i=-1;
        T w=0;
        TV V;
    };
    ARRAY<NODE_DATA> node_data;
    struct INTERACTION_DATA
    {
        int flat_index;
        int node_data_index;
        T w;
    };
    ARRAY<ARRAY<ARRAY<INTERACTION_DATA> > > interaction_data(example.collision_objects_reflection.m);
    for(int i=0;i<example.collision_objects_reflection.m;i++)
    {
        auto* co=example.collision_objects(i);
        auto* cor=example.collision_objects_reflection(i);
        FRAME<TV> frame;
        if(co->func_frame) frame=co->func_frame(example.time+example.dt);
        if(!cor)
        {
            cor=new typename MPM_EXAMPLE<TV>::REFLECT_OBJECT_DATA;
            example.collision_objects_reflection(i)=cor;
            Sample_Reflection_Collision_Object(i);
        }
        interaction_data(i).Resize(cor->X.m);
        typename PARTICLE_GRID_ITERATOR<TV>::SCRATCH scratch;
        for(int p=0;p<cor->X.m;p++)
        {
            auto& id=interaction_data(i)(p);
            TV X=frame*cor->X(p);
            for(PARTICLE_GRID_ITERATOR<TV> it(example.weights,X,false,scratch);it.Valid();it.Next())
            {
                int index=example.mass.Standard_Index(it.Index());
                T w=it.Weight();
                if(!example.mass.array(index)) continue;
                if(w<std::numeric_limits<T>::epsilon()) continue;
                LOG::printf("w %P\n",w);
                int& cn=example.colliding_nodes.array(index);
                if(cn<0)
                {
                    cn=node_data.Append({index,w,example.velocity.array(index)});
                    example.velocity.array(index)=TV();
                }
                else node_data(cn).w+=w;
                id.Append({index,cn,w});
                Add_Debug_Object(VECTOR<TV,2>(X,example.grid.Center(it.Index())),VECTOR<T,3>(1,0,0));
            }
            Add_Debug_Particle(X,VECTOR<T,3>(1,!id.m,0));
        }
    }
    for(int i=0;i<example.collision_objects_reflection.m;i++)
    {
        auto* co=example.collision_objects(i);
        auto* cor=example.collision_objects_reflection(i);
        const auto& id_o=interaction_data(i);
        FRAME<TV> frame;
        if(co->func_frame) frame=co->func_frame(example.time+example.dt);
        for(int p=0;p<id_o.m;p++)
        {
            const auto& id_p=id_o(p);
            if(!id_p.m) continue;
            T mp=0,beta=0;
            TV vp;
            for(const auto& id_i:id_p)
            {
                const NODE_DATA& d=node_data(id_i.node_data_index);
                T w=id_i.w,a=w/d.w,mi=example.mass.array(id_i.flat_index),a_mi=a*mi,w_a_mi=w*a_mi;
                mp+=w_a_mi;
                beta+=w*w_a_mi;
                vp+=w_a_mi*d.V;
            }
            vp/=mp;
            beta/=mp;
            TV X=frame*cor->X(p);
            TV n=co->Normal(X,example.time+example.dt);
            TV b=co->Velocity(X,example.time+example.dt);
            TV vp_b=(vp-b)/-beta;
            T vp_b_n=0;
            if(co->type!=co->stick)
            {
                LOG::printf("FAIL!!!\n");
                vp_b_n=vp_b.Dot(n);
                vp_b=vp_b_n*n;
            }
            TV vp2;
            for(const auto& id_i:id_p)
            {
                const NODE_DATA& d=node_data(id_i.node_data_index);
                T w=id_i.w,a=w/d.w,mi=example.mass.array(id_i.flat_index),a_mi=a*mi,w_a_mi=w*a_mi;
                TV j_ip=w_a_mi*vp_b;
                TV v0=d.V,v1=v0+w*vp_b;
                if(co->type==co->separate)
                {
                    if(vp_b_n<0) v1=v0;
                    else
                    {
                        TV vt=(b-v1).Projected_Orthogonal_To_Unit_Direction(n);
                        T vt_m=vt.Normalize();
                        v1+=min(co->friction*w*vp_b_n,vt_m)*vt;
                    }
                }
                vp2+=w_a_mi*v1;
                LOG::printf("vp2 %P %P %P\n",w,w_a_mi,v1);
                example.velocity.array(id_i.flat_index)+=a*v1;
            }
            vp2/=mp;
            LOG::printf("new V %P\n",vp2);
        }
    }
    for(const auto& cn:node_data)
        example.colliding_nodes.array(cn.i)=-1;
}
//#####################################################################
// Function Sample_Reflection_Collision_Object
//#####################################################################
template<class TV> void MPM_DRIVER<TV>::
Sample_Reflection_Collision_Object(int i)
{
    POISSON_DISK_SURFACE<TV> pds;
    pds.Set_Distance((T).5*example.grid.dX.Min());
    pds.Sample(example.random,example.collision_objects(i)->io,example.collision_objects_reflection(i)->X);
}
//#####################################################################
// Function Step
//#####################################################################
template<class TV> void MPM_DRIVER<TV>::
Step(std::function<void()> func,const char* name,bool dump_substep,bool do_step)
{
    TIMER_SCOPE_FUNC;
    auto& p=example.time_step_callbacks.Get_Or_Insert(name);
    p.x=true; // Flag the callback name as recognized, for sanity checking later
    if(!do_step) return;
    for(int i=0;i<p.y(1).m;i++) p.y(1)(i)();
    func();
    // Note: p may be invalidated by func(), so we must re-access it here.
    const auto& q=example.time_step_callbacks.Get_Or_Insert(name);
    for(int i=0;i<q.y(0).m;i++) q.y(0)(i)();
    if(dump_substep) PHYSBAM_DEBUG_WRITE_SUBSTEP(name,1);
}
//#####################################################################
namespace PhysBAM{
template class MPM_DRIVER<VECTOR<float,2> >;
template class MPM_DRIVER<VECTOR<float,3> >;
template class MPM_DRIVER<VECTOR<double,2> >;
template class MPM_DRIVER<VECTOR<double,3> >;
}
