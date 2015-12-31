//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Grids_Uniform/CELL_ITERATOR.h>
#include <Tools/Grids_Uniform/FACE_ITERATOR.h>
#include <Tools/Krylov_Solvers/CONJUGATE_GRADIENT.h>
#include <Tools/Krylov_Solvers/MINRES.h>
#include <Tools/Log/DEBUG_SUBSTEPS.h>
#include <Tools/Log/LOG.h>
#include <Tools/Log/SCOPE.h>
#include <Tools/Nonlinear_Equations/NEWTONS_METHOD.h>
#include <Deformables/Collisions_And_Interactions/IMPLICIT_OBJECT_COLLISION_PENALTY_FORCES.h>
#include <Deformables/Constitutive_Models/MPM_PLASTICITY_MODEL.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_DRIVER.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_EXAMPLE.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_PARTICLES.h>
#include <Hybrid_Methods/Forces/MPM_PLASTIC_FINITE_ELEMENTS.h>
#include <Hybrid_Methods/Iterators/GATHER_SCATTER.h>
#include <Hybrid_Methods/Iterators/PARTICLE_GRID_WEIGHTS.h>
#include <Hybrid_Methods/System/FLUID_KRYLOV_SYSTEM.h>
#include <Hybrid_Methods/System/FLUID_KRYLOV_VECTOR.h>
#include <Hybrid_Methods/System/KKT_KRYLOV_SYSTEM.h>
#include <Hybrid_Methods/System/KKT_KRYLOV_VECTOR.h>
#include <Hybrid_Methods/System/MPM_KRYLOV_VECTOR.h>
#include <Hybrid_Methods/System/MPM_OBJECTIVE.h>
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
    rhs(*new MPM_KRYLOV_VECTOR<TV>(example.valid_grid_indices)),
    fluid_sys(*new FLUID_KRYLOV_SYSTEM<TV>(example)),
    fluid_p(*new FLUID_KRYLOV_VECTOR<TV>(example.valid_pressure_dofs)),
    fluid_rhs(*new FLUID_KRYLOV_VECTOR<TV>(example.valid_pressure_dofs)),
    kkt_sys(*new KKT_KRYLOV_SYSTEM<TV>(example)),
    kkt_lhs(*new KKT_KRYLOV_VECTOR<TV>(example.valid_pressure_indices)),
    kkt_rhs(*new KKT_KRYLOV_VECTOR<TV>(example.valid_pressure_indices))
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
    delete &fluid_sys;
    delete &fluid_p;
    delete &fluid_rhs;
    av.Delete_Pointers_And_Clean_Memory();
    bv.Delete_Pointers_And_Clean_Memory();
    cv.Delete_Pointers_And_Clean_Memory();
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

    example.particles.Store_B(example.use_affine && !(example.incompressible || example.kkt));
    example.particles.Store_C(example.use_affine && (example.incompressible || example.kkt));
    example.particles.Store_S(example.use_oldroyd);
    example.current_velocity=&example.velocity;
    PHYSBAM_ASSERT(!example.particles.store_B || !example.particles.store_C);

    if(example.incompressible){
        fluid_p.p.Resize(example.grid.Domain_Indices(example.ghost));
        fluid_rhs.p.Resize(example.grid.Domain_Indices(example.ghost));
        example.mass_f.Resize(example.grid.Domain_Indices(example.ghost));
        example.volume_f.Resize(example.grid.Domain_Indices(example.ghost));
        example.density_f.Resize(example.grid.Domain_Indices(example.ghost));
        example.velocity_f.Resize(example.grid.Domain_Indices(example.ghost));
        example.velocity_new_f.Resize(example.grid.Domain_Indices(example.ghost));}
    else if(example.kkt){
        kkt_lhs.u.Resize(example.grid.Domain_Indices(example.ghost));
        kkt_lhs.p.Resize(example.grid.Domain_Indices(example.ghost));
        kkt_rhs.u.Resize(example.grid.Domain_Indices(example.ghost));
        kkt_rhs.p.Resize(example.grid.Domain_Indices(example.ghost));
        example.one_over_lambda.Resize(example.grid.Domain_Indices(example.ghost));
        example.J.Resize(example.grid.Domain_Indices(example.ghost));
        example.density.Resize(example.grid.Domain_Indices(example.ghost));}
    if(example.incompressible || example.kkt){
        example.volume.Resize(example.grid.Domain_Indices(example.ghost));
        example.cell_solid.Resize(example.grid.Domain_Indices(example.ghost));
        example.cell_pressure.Resize(example.grid.Domain_Indices(example.ghost));
        example.cell_C.Resize(example.grid.Domain_Indices(example.ghost));
        for(CELL_ITERATOR<TV> iterator(example.grid,example.ghost,GRID<TV>::GHOST_REGION);iterator.Valid();iterator.Next())
            example.cell_solid(iterator.Cell_Index())=true;
        if(example.kkt){
            RANGE<TV_INT> range(example.grid.Cell_Indices(0));
            TV_INT max_corner=range.max_corner;
            LOG::printf("max_corner=%P\n",max_corner);
            for(CELL_ITERATOR<TV> iterator(example.grid,example.ghost,GRID<TV>::WHOLE_REGION);iterator.Valid();iterator.Next()){
                for(int a=0;a<TV::m;a++){
                    int max_corner_axis=max_corner(a)-TV_INT::Axis_Vector(a)(a);
                    if(iterator.Cell_Index()(a)==0 || iterator.Cell_Index()(a)==1 || iterator.Cell_Index()(a)==max_corner_axis || iterator.Cell_Index()(a)==max_corner_axis || iterator.Cell_Index()(a)==max_corner_axis-1)
                        example.cell_solid(iterator.Cell_Index())=true;}}}}

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

    Update_Simulated_Particles();
    Print_Particle_Stats("particle state",example.dt);
    Update_Particle_Weights();
    example.gather_scatter.Prepare_Scatter(example.particles);
    Particle_To_Grid();
    Print_Grid_Stats("after particle to grid",example.dt,example.velocity,0);
    Print_Energy_Stats("after particle to grid",example.velocity);
    PHYSBAM_DEBUG_WRITE_SUBSTEP("after particle to grid",0,1);
    Apply_Forces();
    Print_Grid_Stats("after forces",example.dt,example.velocity_new,&example.velocity);
    PHYSBAM_DEBUG_WRITE_SUBSTEP("after forces",0,1);
    if(example.kkt)
        Solve_KKT_System();
    else if(example.incompressible)
        Make_Incompressible();
    else
        Grid_To_Particle();
    PHYSBAM_DEBUG_WRITE_SUBSTEP("after grid to particle",0,1);
    Update_Plasticity_And_Hardening();

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
            T next_time=example.time+example.dt;
            if(next_time>time_at_frame){
                next_time=time_at_frame;
                done=true;}
            else if(next_time+example.dt>time_at_frame) next_time=(example.time+time_at_frame)/2;
            example.dt=next_time-example.time;
            LOG::cout<<"substep dt: "<<example.dt<<std::endl;

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
    LOG::SCOPE scope("Write_Output_Files");
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
    if(example.incompressible)
        for(int i=0;i<TV::m;++i)
            example.face_weights(i)->Update(example.particles.X);
}
//#####################################################################
// Function Particle_To_Grid
//#####################################################################
template<class TV> void MPM_DRIVER<TV>::
Particle_To_Grid()
{
    MPM_PARTICLES<TV>& particles=example.particles;
    example.current_velocity=&example.velocity;

#pragma omp parallel for
    for(int i=0;i<example.mass.array.m;i++){
        example.mass.array(i)=0;
        if(example.incompressible || example.kkt) example.volume.array(i)=0;
        example.velocity.array(i)=TV();
        example.velocity_new.array(i)=TV();
        if(example.kkt){
            example.one_over_lambda.array(i)=0;
            example.J.array(i)=0;}}

    T scale=1;
    if(particles.store_B && !example.weights->use_gradient_transfer)
        scale=example.weights->Constant_Scalar_Inverse_Dp();
    bool use_gradient=!example.incompressible && !example.kkt && example.weights->use_gradient_transfer;
    ARRAY_VIEW<MATRIX<T,TV::m> > dV((example.incompressible || example.kkt)?particles.C:particles.B);
    example.gather_scatter.template Scatter<int>(true,0,
        [this,scale,&particles,use_gradient,dV](int p,const PARTICLE_GRID_ITERATOR<TV>& it,int data)
        {
            T w=it.Weight();
            TV_INT index=it.Index();
            example.mass(index)+=w*particles.mass(p);
            if(example.incompressible || example.kkt) example.volume(index)+=w*particles.volume(p);
            TV V=w*particles.V(p);
            if(example.use_affine){
                if(use_gradient) V+=particles.B(p)*it.Gradient();
                else V+=dV(p)*(w*scale*(example.grid.Center(index)-particles.X(p)));}
            example.velocity(index)+=particles.mass(p)*V;
            if(example.kkt){
                if(particles.lambda(p)!=FLT_MAX) example.one_over_lambda(index)+=w*particles.mass(p)/particles.lambda(p);
                example.J(index)+=w*particles.mass(p)*particles.F(p).Determinant();}
        });

    example.valid_grid_indices.Remove_All();
    example.valid_grid_cell_indices.Remove_All();
    if(example.incompressible || example.kkt){
        example.cell_pressure.Fill(0);
        example.valid_pressure_indices.Remove_All();
        example.valid_pressure_cell_indices.Remove_All();
        example.valid_pressure_dofs.Remove_All();
        example.valid_pressure_cell_dofs.Remove_All();}

    for(RANGE_ITERATOR<TV::m> it(example.mass.domain);it.Valid();it.Next()){
        int i=example.mass.Standard_Index(it.index);
        if(example.mass.array(i)){
            example.valid_grid_indices.Append(i);
            example.valid_grid_cell_indices.Append(it.index);
            example.velocity.array(i)/=example.mass.array(i);
            if(example.incompressible){
                example.cell_pressure.array(i)=1;
                if(!example.cell_solid.array(i)){
                    example.valid_pressure_indices.Append(i);
                    example.valid_pressure_cell_indices.Append(it.index);}}
            if(example.kkt){
                example.one_over_lambda.array(i)/=example.mass.array(i);
                example.J.array(i)/=example.mass.array(i);
                example.density.array(i)=example.mass.array(i)*example.grid.One_Over_Cell_Size();
                if(!example.cell_solid.array(i)){
                    example.cell_pressure.array(i)=1;
                    example.valid_pressure_indices.Append(i);
                    example.valid_pressure_cell_indices.Append(it.index);}
                else{ // zero out velocity, 1/lambda, J in solid region
                    example.velocity.array(i)=TV();
                    example.one_over_lambda.array(i)=(T)0;
                    example.J.array(i)=(T)FLT_MAX;
                    example.density.array(i)=(T)0;}}}
        else example.velocity.array(i)=TV();}
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
        MATRIX<T,TV::m> B,grad_Vp,D;

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
            MATRIX<T,TV::m> A=dt*grad_Vp+1;
            if(example.quad_F_coeff) A+=sqr(dt*grad_Vp)*example.quad_F_coeff;
            particles.F(p)=A*particles.F(p);
            if(particles.store_S){
                T k=example.dt*example.inv_Wi;
                particles.S(p)=(SYMMETRIC_MATRIX<T,TV::m>::Conjugate(A,particles.S(p))+k)/(1+k);}
            if(example.use_affine && example.use_early_gradient_transfer)
                B=grad_Vp/example.weights->Constant_Scalar_Inverse_Dp();
            else if(example.use_affine)
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
                        +MATRIX<T,TV::m>::Outer_Product(Z-particles.X(p)-xi_new+xp_new,V_grid));
                    if(particles.store_C && example.weights->Order()>1) D+=it.Weight()*MATRIX<T,TV::m>::Outer_Product(Z-particles.X(p),Z-particles.X(p));}
            // apply penalty friction on V_pic
            if(example.penalty_map.m){
                const ARRAY<VECTOR<int,2> >& cc=example.penalty_map(p);
                for(int ob=0;ob<cc.m;ob++){
                    IMPLICIT_OBJECT_COLLISION_PENALTY_FORCES<TV>* cf=dynamic_cast<IMPLICIT_OBJECT_COLLISION_PENALTY_FORCES<TV>* >(example.lagrangian_forces(cc(ob).x));
                    int i=cc(ob).y;
                    TV f=cf->grad_pe(i);
                    T fn=f.Normalize();
                    TV t=V_pic.Projected_Orthogonal_To_Unit_Direction(f);
                    T t_mag=t.Normalize();
                    if(t_mag<=cf->coefficient_of_friction*fn/particles.mass(p))
                        V_pic.Project_On_Unit_Direction(f);
                    else V_pic-=cf->coefficient_of_friction/particles.mass(p)*fn*t*example.dt;}}

            particles.V(p)=V_pic;
            Perform_Particle_Collision(p,example.time+example.dt);
            if(particles.store_B) particles.B(p)=B;
            if(particles.store_C) particles.C(p)=example.weights->Order()==1?grad_Vp:B*D.Inverse();
            if(example.use_midpoint) particles.X(p)+=(particles.V(p)+Vn_interpolate)*(dt/2);
            else particles.X(p)+=particles.V(p)*dt;
            particles.V(p)=V_flip*example.flip+V_pic*(1-example.flip);

            if(!example.grid.domain.Lazy_Inside(particles.X(p))) particles.valid(p)=false;}}
}
//#####################################################################
// Function Face_To_Particle
//#####################################################################
template<class TV> void MPM_DRIVER<TV>::
Face_To_Particle()
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
        MATRIX<T,TV::m> C;
            for(int i=0;i<TV::m;++i)
                for(PARTICLE_GRID_ITERATOR<TV> it(example.face_weights(i),p,true,scratch);it.Valid();it.Next()){
            T w=it.Weight();
            FACE_INDEX<TV::m> face_index(i,it.Index());
            T V_grid=example.velocity_new_f(face_index);
            V_pic(i)+=w*V_grid;
            V_flip(i)+=w*(V_grid-example.velocity_f(face_index));
            if(example.use_midpoint){
                PHYSBAM_NOT_IMPLEMENTED("Midpoint with Face_To_Particle is not supported");}
            if(example.use_affine){
                TV tmp=example.velocity_new_f(face_index)*it.Gradient();
                for(int k=0;k<TV::m;++k)
                            C(i,k)+=tmp(k);}}
            particles.V(p)=V_pic;
            Perform_Particle_Collision(p,example.time+example.dt);
            if(example.use_midpoint) PHYSBAM_NOT_IMPLEMENTED("Midpoint with Face_To_Particle is not supported");
            else particles.X(p)+=particles.V(p)*dt;
            particles.V(p)=V_flip*example.flip+V_pic*(1-example.flip);
            if(particles.store_C) particles.C(p)=C;
            particles.F(p)+=dt*C*particles.F(p);
            if(!example.grid.domain.Lazy_Inside(particles.X(p))) particles.valid(p)=false;}}
}
//#####################################################################
// Function Compute_Cell_C
//#####################################################################
template<class TV> void MPM_DRIVER<TV>::
Compute_Cell_C()
{
    example.cell_C.Fill(MATRIX<T,TV::m>());
    T one_over_dx=example.grid.one_over_dX(0);
#pragma omp parallel for
    for(int t=0;t<example.valid_grid_cell_indices.m;t++){
        TV_INT index=example.valid_grid_cell_indices(t);
        MATRIX<T,TV::m> C;
        for(int i=0;i<TV::m;i++)
            for(int j=0;j<TV::m;j++){
                TV_INT e=TV_INT::Axis_Vector(j);
                if(example.mass(index+e)) 
                    C(i,j)=(example.velocity_new(index+e)(i)-example.velocity_new(index)(i))*one_over_dx;
                else if(example.mass(index-e)) 
                    C(i,j)=(example.velocity_new(index)(i)-example.velocity_new(index-e)(i))*one_over_dx;
                else
                    PHYSBAM_FATAL_ERROR("Failed to compute cell centered C");}
        example.cell_C(index)=C;}
}
//#####################################################################
// Function Cell_To_Face
//#####################################################################
template<class TV> void MPM_DRIVER<TV>::
Cell_To_Face()
{
    // TODO: parallelize
    example.mass_f.Fill(0);
    example.velocity_f.Fill(0);
    example.velocity_new_f.Fill(0);

#pragma omp parallel for
    for(int t=0;t<example.valid_grid_cell_indices.m;t++){
        TV_INT index=example.valid_grid_cell_indices(t);
        for(int a=0;a<TV::m;a++){
            FACE_INDEX<TV::m> face(a,example.grid.First_Face_Index_In_Cell(a,index));
            example.mass_f(face)+=example.mass(index);
            example.velocity_f(face)+=example.mass(index)*example.velocity_new(index)(a);}}
#pragma omp parallel for
    for(int t=0;t<example.valid_grid_cell_indices.m;t++){
        TV_INT index=example.valid_grid_cell_indices(t);
        for(int a=0;a<TV::m;a++){
            FACE_INDEX<TV::m> face(a,example.grid.Second_Face_Index_In_Cell(a,index));
            example.mass_f(face)+=example.mass(index);
            example.velocity_f(face)+=example.mass(index)*example.velocity_new(index)(a);}}

    // TODO: improve efficiency
    for(FACE_ITERATOR<TV> iterator(example.grid,example.ghost);iterator.Valid();iterator.Next()){
        FACE_INDEX<TV::m> face=iterator.Full_Index();
        if(example.mass_f(face)){
            example.velocity_f(face)/=example.mass_f(face);
            example.mass_f(face)*=(T).5;}}
}
//#####################################################################
// Function Cell_To_Face_C
//#####################################################################
template<class TV> void MPM_DRIVER<TV>::
Cell_To_Face_C()
{
    // TODO: parallelize
    example.mass_f.Fill(0);
    example.volume_f.Fill(0);
    example.velocity_f.Fill(0);
    example.velocity_new_f.Fill(0);
    example.density_f.Fill(0);

#pragma omp parallel for
    for(int t=0;t<example.valid_grid_cell_indices.m;t++){
        TV_INT index=example.valid_grid_cell_indices(t);
        for(int a=0;a<TV::m;a++){
            FACE_INDEX<TV::m> face(a,example.grid.First_Face_Index_In_Cell(a,index));
            example.mass_f(face)+=example.mass(index);
            example.volume_f(face)+=example.volume(index);
            example.velocity_f(face)+=example.mass(index)*(example.velocity_new(index)(a)-(example.cell_C(index)*TV::Axis_Vector(a)*example.grid.dX(a))(a)*(T).5);}}
#pragma omp parallel for
    for(int t=0;t<example.valid_grid_cell_indices.m;t++){
        TV_INT index=example.valid_grid_cell_indices(t);
        for(int a=0;a<TV::m;a++){
            FACE_INDEX<TV::m> face(a,example.grid.Second_Face_Index_In_Cell(a,index));
            example.mass_f(face)+=example.mass(index);
            example.volume_f(face)+=example.volume(index);
            example.velocity_f(face)+=example.mass(index)*(example.velocity_new(index)(a)+(example.cell_C(index)*TV::Axis_Vector(a)*example.grid.dX(a))(a)*(T).5);}}

    // TODO: improve efficiency
    for(FACE_ITERATOR<TV> iterator(example.grid,example.ghost);iterator.Valid();iterator.Next()){
        FACE_INDEX<TV::m> face=iterator.Full_Index();
        if(example.mass_f(face)){
            example.velocity_f(face)/=example.mass_f(face);
            example.mass_f(face)*=(T).5;
            example.volume_f(face)*=(T).5;
            example.density_f(face)=example.mass_f(face)/example.volume_f(face);}}
}
//#####################################################################
// Function Face_To_Cell
//#####################################################################
template<class TV> void MPM_DRIVER<TV>::
Face_To_Cell()
{
#pragma omp parallel for
    for(int t=0;t<example.valid_grid_cell_indices.m;t++){
        TV_INT index=example.valid_grid_cell_indices(t);
        for(int a=0;a<TV::m;a++){
            FACE_INDEX<TV::m> faceF(a,example.grid.First_Face_Index_In_Cell(a,index));
            FACE_INDEX<TV::m> faceS(a,example.grid.Second_Face_Index_In_Cell(a,index));
            example.velocity_new(index)(a)=(example.velocity_new_f(faceF)+example.velocity_new_f(faceS))*(T).5;}}
}
//#####################################################################
// Function Make_Incompressible
//#####################################################################
template<class TV> void MPM_DRIVER<TV>::
Make_Incompressible()
{
    Compute_Cell_C();
    Cell_To_Face_C();
    example.velocity_new_f=example.velocity_f;
    Pressure_Projection();
    if(example.use_f2p)
        Face_To_Particle();
    else{
        Face_To_Cell();
        Grid_To_Particle();}
}
//#####################################################################
// Function Solve_KKT_System
//#####################################################################
template<class TV> void MPM_DRIVER<TV>::
Solve_KKT_System()
{
    // Zero out variables in solid
    for(int i=0;i<example.valid_grid_indices.m;i++){
        TV_INT index=example.valid_grid_cell_indices(i);
        if(example.cell_solid(index)){
            example.mass(index)=0;
            example.density(index)=0;
            example.velocity_new(index)=TV();}}
    // Construct DT
    T one_over_dt=(T)1/example.dt;
    int d=TV::m;
    ARRAY<T,TV_INT>& m=example.mass;
    example.DT.Reset(d*example.velocity.array.m);
    for(int row=0;row<example.valid_pressure_indices.m;row++){
        TV_INT id=example.valid_pressure_cell_indices(row);
        for(int a=0;a<d;a++){
            TV_INT ej=TV_INT::Axis_Vector(a);
            for(int i=-1;i<3;i++) example.DT.Append_Entry_To_Current_Row(d*example.velocity.Standard_Index(id+i*ej)+a,(T)0);}
        example.DT.Finish_Row();}
    example.DT.Sort_Entries();
    for(int row=0;row<example.valid_pressure_indices.m;row++){
        TV_INT id=example.valid_pressure_cell_indices(row);
        for(int a=0;a<TV::m;a++){
            TV_INT ej=TV_INT::Axis_Vector(a);
            // Go through face A and face B
            for(int f=0;f<2;f++){
                int sign=2*f-1;
                TV_INT si=id+sign*ej;
                if(!example.cell_solid(si)){
                    example.DT(row,d*example.velocity.Standard_Index(id)+a)+=sign*m(id)/(m(id)+m(si));
                    example.DT(row,d*example.velocity.Standard_Index(si)+a)+=sign*m(si)/(m(id)+m(si));
                    Add_C_Contribution_To_DT(id,a,row,m(id)/(2*(m(id)+m(si))));
                    Add_C_Contribution_To_DT(si,a,row,-m(si)/(2*(m(id)+m(si))));}}}}
    example.DT*=example.grid.one_over_dX(0);

    // Construct lhs and rhs
    kkt_lhs.u.Fill(TV());kkt_rhs.u.Fill(TV());
    kkt_lhs.p.Fill((T)0);kkt_rhs.p.Fill((T)0);
    for(int t=0;t<example.valid_pressure_cell_indices.m;t++){
        TV_INT valid_cell=example.valid_pressure_cell_indices(t);
        kkt_rhs.u(valid_cell)=example.velocity_new(valid_cell)*sqrt(example.density(valid_cell)/example.dt);
        kkt_rhs.p(valid_cell)=one_over_dt*((T)1-(T)1/example.J(valid_cell));}
    // MINRES 
    MINRES<T> mr;
    int max_iterations=100000;
    mr.Solve(kkt_sys,kkt_lhs,kkt_rhs,cv,1e-12,0,max_iterations);
    // Update velocity on grid
    example.velocity_new=kkt_lhs.u;
    for(int i=0;i<example.valid_pressure_indices.m;i++){
        TV_INT index=example.valid_pressure_cell_indices(i);
        example.velocity_new(index)*=sqrt(example.dt/example.density(index));}
    // Apply BC
    for(int t=0; t<example.valid_grid_indices.m; t++){
        TV_INT index=example.valid_grid_cell_indices(t);
        if(example.cell_solid(index)) for(int b=0;b<TV::m;b++) for(int a=0;a<TV::m; a++){
            const TV_INT axis=TV_INT::Axis_Vector(a);
            if(!example.cell_solid(index+axis)&&!example.cell_solid(index+2*axis))
                example.velocity_new(index)(b)=a==b?-example.velocity_new(index+axis)(b):example.velocity_new(index+axis)(b);
            else if(example.cell_solid(index+axis)&&!example.cell_solid(index+2*axis))
                example.velocity_new(index)(b)=a==b?-example.velocity_new(index+2*axis)(b):example.velocity_new(index+2*axis)(b);
            else if(!example.cell_solid(index-axis)&&!example.cell_solid(index-2*axis))
                example.velocity_new(index)(b)=a==b?-example.velocity_new(index-axis)(b):example.velocity_new(index-axis)(b);
            else if(example.cell_solid(index-axis)&&!example.cell_solid(index-2*axis))
                example.velocity_new(index)(b)=a==b?-example.velocity_new(index-2*axis)(b):example.velocity_new(index-2*axis)(b);}}
    // Transfer to particles
    Grid_To_Particle();
    PHYSBAM_DEBUG_WRITE_SUBSTEP("after minres in kkt",0,1);
}
//#####################################################################
// Function Update_Plasticity_And_Hardening
//#####################################################################
template<class TV> void MPM_DRIVER<TV>::
Update_Plasticity_And_Hardening()
{
    for(int i=0;i<example.plasticity_models.m;i++)
        example.plasticity_models(i)->Update_Particles();
}
//#####################################################################
// Function Add_C_Contribution_To_DT
//#####################################################################
template<class TV> void MPM_DRIVER<TV>::
Add_C_Contribution_To_DT(TV_INT id,int a,int row,T val)
{
    const TV_INT axis=TV_INT::Axis_Vector(a);
    TV_INT right_index=example.cell_pressure(id+axis)?id+axis:id;
    int d=TV::m;
    example.DT(row,d*example.velocity.Standard_Index(right_index)+a)+=val;
    example.DT(row,d*example.velocity.Standard_Index(right_index-axis)+a)-=val;
}
//#####################################################################
// Function Pressure_Projection
//#####################################################################
template<class TV> typename TV::SCALAR MPM_DRIVER<TV>::
Pressure_Projection()
{
    CONJUGATE_GRADIENT<T>* cg=new CONJUGATE_GRADIENT<T>();
    int max_solver_iterations=1000;
    cg->finish_before_indefiniteness=true;
    cg->relative_tolerance=false;
    T one_over_dx=example.grid.one_over_dX(0);
    
    // update valid dofs
    fluid_p.p.Fill(0);
    fluid_rhs.p.Fill(0);
    ARRAY<TV_INT> boundary_correction,dofs_near_boundary;
    boundary_correction.Remove_All(); dofs_near_boundary.Remove_All();
    auto active_dofs=[&](const TV_INT& index){
        bool not_near=true;
        for(int a=0;a<TV_INT::m;++a){
            const TV_INT axis=TV_INT::Axis_Vector(a);
            not_near&=!(example.cell_solid(index-axis)||example.cell_solid(index+axis)||example.cell_solid(index-2*axis)||example.cell_solid(index+2*axis));}
        return not_near&&example.cell_pressure(index);};
    for(int t=0;t<example.valid_pressure_cell_indices.m;t++){
        TV_INT index=example.valid_pressure_cell_indices(t);
        if(active_dofs(index)){
            example.valid_pressure_dofs.Append(example.valid_pressure_indices(t));
            example.valid_pressure_cell_dofs.Append(index);}
        else
            boundary_correction.Append(index);
        for(int a=0;a<TV::m;a++){
            const TV_INT axis=TV_INT::Axis_Vector(a);
            if((example.cell_solid(index-3*axis)&&!example.cell_solid(index-2*axis))||(example.cell_solid(index+3*axis)&&!example.cell_solid(index+2*axis)))
                dofs_near_boundary.Append(index);}}
    for(RANGE_ITERATOR<TV::m> it(example.mass.domain);it.Valid();it.Next()){
        int i=example.mass.Standard_Index(it.index);
        if(example.mass.array(i)&&active_dofs(it.index)){
            example.cell_pressure.array(i)=2;}}
    for(int t=0;t<example.valid_pressure_cell_dofs.m;t++){
        TV_INT index=example.valid_pressure_cell_dofs(t);
        for(int a=0;a<TV::m;a++){
            const TV_INT axis=TV_INT::Axis_Vector(a);
            FACE_INDEX<TV::m> faceF(a,example.grid.First_Face_Index_In_Cell(a,index));
            FACE_INDEX<TV::m> faceS(a,example.grid.Second_Face_Index_In_Cell(a,index));
            fluid_rhs.p(index)+=example.velocity_f(faceF)-example.velocity_f(faceS);
            if(example.cell_solid(index-3*axis)&&!example.cell_solid(index-2*axis))
                fluid_rhs.p(index)-=example.velocity_f(faceF);
            if(example.cell_solid(index+3*axis)&&!example.cell_solid(index+2*axis))
                fluid_rhs.p(index)+=example.velocity_f(faceS);}
        fluid_rhs.p(index)*=one_over_dx;}

    bool converged=cg->Solve(fluid_sys,fluid_p,fluid_rhs,bv,1e-8,0,max_solver_iterations);
    if(!converged) LOG::printf("SOLVER DID NOT CONVERGE.\n");
    example.velocity_new_f=example.velocity_f;
    
    // update velocity on faces
    for(int t=0;t<example.valid_pressure_cell_dofs.m;t++){
        TV_INT index=example.valid_pressure_cell_dofs(t);
        for(int a=0;a<TV::m;a++){
            FACE_INDEX<TV::m> faceF(a,example.grid.First_Face_Index_In_Cell(a,index));
            FACE_INDEX<TV::m> faceS(a,example.grid.Second_Face_Index_In_Cell(a,index));
            example.velocity_new_f(faceF)-=fluid_p.p(index)*one_over_dx/(T)example.density_f(faceF);
            example.velocity_new_f(faceS)+=fluid_p.p(index)*one_over_dx/(T)example.density_f(faceS);}}

    for(int t=0;t<dofs_near_boundary.m;t++){
        TV_INT index=dofs_near_boundary(t);
        for(int a=0;a<TV::m;a++){
            const TV_INT axis=TV_INT::Axis_Vector(a);
            if(example.cell_solid(index-3*axis)&&!example.cell_solid(index-2*axis)){
                FACE_INDEX<TV::m> faceF(a,example.grid.First_Face_Index_In_Cell(a,index));
                example.velocity_new_f(faceF)=0;}
            if(example.cell_solid(index+3*axis)&&!example.cell_solid(index+2*axis)){
                FACE_INDEX<TV::m> faceS(a,example.grid.Second_Face_Index_In_Cell(a,index));
                example.velocity_new_f(faceS)=0;}}}
            
    for(int t=0;t<boundary_correction.m;t++){
        TV_INT index=boundary_correction(t);
        for(int a=0;a<TV::m;a++){
            const TV_INT axis=TV_INT::Axis_Vector(a);
            if(example.cell_solid(index-2*axis)&&!example.cell_solid(index-axis)){
                for(int b=0;b<TV::m;b++){
                    FACE_INDEX<TV::m> faceF(b,example.grid.First_Face_Index_In_Cell(b,index));
                    FACE_INDEX<TV::m> faceS(b,example.grid.Second_Face_Index_In_Cell(b,index));
                    if(b==a){
                        example.velocity_new_f(faceF)=-example.velocity_new_f(FACE_INDEX<TV::m>(b,example.grid.Second_Face_Index_In_Cell(b,index+axis)));}
                    else{
                        example.velocity_new_f(faceF)=example.velocity_new_f(FACE_INDEX<TV::m>(b,example.grid.First_Face_Index_In_Cell(b,index+axis)));
                        example.velocity_new_f(faceS)=example.velocity_new_f(FACE_INDEX<TV::m>(b,example.grid.Second_Face_Index_In_Cell(b,index+axis)));}}}
            if(example.cell_solid(index+2*axis)&&!example.cell_solid(index+axis)){
                for(int b=0;b<TV::m;b++){
                    FACE_INDEX<TV::m> faceF(b,example.grid.First_Face_Index_In_Cell(b,index));
                    FACE_INDEX<TV::m> faceS(b,example.grid.Second_Face_Index_In_Cell(b,index));
                    if(b==a){
                        example.velocity_new_f(faceS)=-example.velocity_new_f(FACE_INDEX<TV::m>(b,example.grid.First_Face_Index_In_Cell(b,index-axis)));}
                    else{
                        example.velocity_new_f(faceF)=example.velocity_new_f(FACE_INDEX<TV::m>(b,example.grid.First_Face_Index_In_Cell(b,index-axis)));
                        example.velocity_new_f(faceS)=example.velocity_new_f(FACE_INDEX<TV::m>(b,example.grid.Second_Face_Index_In_Cell(b,index-axis)));}}}}}

    for(int t=0;t<boundary_correction.m;t++){
        TV_INT index=boundary_correction(t);
        for(int a=0;a<TV::m;a++){
            const TV_INT axis=TV_INT::Axis_Vector(a);
            if(example.cell_solid(index-2*axis)&&!example.cell_solid(index-axis)){
                for(int b=0;b<TV::m;b++){
                    FACE_INDEX<TV::m> faceF(b,example.grid.First_Face_Index_In_Cell(b,index));
                    FACE_INDEX<TV::m> faceS(b,example.grid.Second_Face_Index_In_Cell(b,index));
                    if(b==a){
                        example.velocity_new_f(faceF)=-example.velocity_new_f(FACE_INDEX<TV::m>(b,example.grid.Second_Face_Index_In_Cell(b,index+axis)));}
                    else{
                        example.velocity_new_f(faceF)=example.velocity_new_f(FACE_INDEX<TV::m>(b,example.grid.First_Face_Index_In_Cell(b,index+axis)));
                        example.velocity_new_f(faceS)=example.velocity_new_f(FACE_INDEX<TV::m>(b,example.grid.Second_Face_Index_In_Cell(b,index+axis)));}}}
            if(example.cell_solid(index+2*axis)&&!example.cell_solid(index+axis)){
                for(int b=0;b<TV::m;b++){
                    FACE_INDEX<TV::m> faceF(b,example.grid.First_Face_Index_In_Cell(b,index));
                    FACE_INDEX<TV::m> faceS(b,example.grid.Second_Face_Index_In_Cell(b,index));
                    if(b==a){
                        example.velocity_new_f(faceS)=-example.velocity_new_f(FACE_INDEX<TV::m>(b,example.grid.First_Face_Index_In_Cell(b,index-axis)));}
                    else{
                        example.velocity_new_f(faceF)=example.velocity_new_f(FACE_INDEX<TV::m>(b,example.grid.First_Face_Index_In_Cell(b,index-axis)));
                        example.velocity_new_f(faceS)=example.velocity_new_f(FACE_INDEX<TV::m>(b,example.grid.Second_Face_Index_In_Cell(b,index-axis)));}}}}}

    // Diagnostics
    T divmax=0;
    TV_INT where;
    int cell_type;
    for(int t=0;t<example.valid_pressure_cell_indices.m;t++){
        TV_INT index=example.valid_pressure_cell_indices(t);
        T div=0;
        for(int a=0;a<TV::m;a++){
            const TV_INT axis=TV_INT::Axis_Vector(a);
            FACE_INDEX<TV::m> faceF(a,example.grid.First_Face_Index_In_Cell(a,index));
            FACE_INDEX<TV::m> faceS(a,example.grid.Second_Face_Index_In_Cell(a,index));
            div+=example.velocity_new_f(faceS)-example.velocity_new_f(faceF);}
        if(abs(div)>divmax) {divmax=abs(div); where=index; cell_type=example.cell_pressure(index);}}
    divmax*=one_over_dx;
    LOG::cout<<"Divergence after projection\t"<<divmax<<std::endl;
    LOG::printf("where: %P\t cell_type: %d\n",where,cell_type);
    return divmax;
}
//#####################################################################
// Function Apply_Forces
//#####################################################################
template<class TV> void MPM_DRIVER<TV>::
Apply_Forces()
{
    example.Capture_Stress();
    objective.Reset();
    LOG::printf("max velocity: %P\n",Max_Particle_Speed());
    if(example.use_symplectic_euler){
        objective.tmp2*=0;
        example.Precompute_Forces(example.time,example.dt,false);
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
        if(example.asymmetric_system){
            newtons_method.use_gmres=true;
            newtons_method.use_cg=false;
            newtons_method.Make_Vanilla_Newton();}

        example.Update_Lagged_Forces(example.time);
        newtons_method.require_one_iteration=!objective.Initial_Guess(dv,newtons_method.tolerance,example.asymmetric_system);
        if(example.test_diff) objective.Test_Diff(dv);

        objective.system.forced_collisions.Remove_All();
        bool converged=newtons_method.Newtons_Method(objective,objective.system,dv,av);
        if(!converged) LOG::cout<<"WARNING: Newton's method did not converge"<<std::endl;}

    Apply_Friction();
    objective.Restore_F();

#pragma omp parallel for
    for(int i=0;i<example.valid_grid_indices.m;i++){
        int j=example.valid_grid_indices(i);
        example.velocity_new.array(j)=dv.u.array(j)+objective.v0.u.array(j);}
    example.velocity_new.array.Subset(objective.system.stuck_nodes)=objective.system.stuck_velocity;
    example.current_velocity=&example.velocity_new;

    Register_Active_Penalty_Collisions();
}
//#####################################################################
// Function Register_Active_Penalty_Collisions
//#####################################################################
template<class TV> void MPM_DRIVER<TV>::
Register_Active_Penalty_Collisions()
{
    ARRAY<ARRAY<VECTOR<int,2> > >& map=example.penalty_map;
    MPM_PARTICLES<TV>& particles=example.particles;
    map.Clean_Memory();
    map.Resize(particles.number);
    for(int ii=0;ii<example.lagrangian_forces.m;ii++)
        if(IMPLICIT_OBJECT_COLLISION_PENALTY_FORCES<TV>* cf=dynamic_cast<IMPLICIT_OBJECT_COLLISION_PENALTY_FORCES<TV>* >(example.lagrangian_forces(ii)))
            for(int i=0;i<cf->penetrating_particles.m;i++){
                int p=cf->penetrating_particles(i);
                map(p).Append(VECTOR<int,2>(ii,i));}
}
//#####################################################################
// Function Perform_Particle_Collision
//#####################################################################
template<class TV> void MPM_DRIVER<TV>::
Perform_Particle_Collision(int p,T time)
{
    if(!example.use_particle_collision) return;;
    for(int i=0;i<example.collision_objects.m;i++){
        TV X=example.particles.X(p);
        MPM_COLLISION_OBJECT<TV>* io=example.collision_objects(i);
        T phi=io->Phi(X,time);
        if(phi>=0) continue;
        if(example.collision_objects(i)->type==COLLISION_TYPE::stick) return;
        X-=phi*io->Normal(X,time);
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
        T coefficient_of_friction=example.collision_objects(c.object)->friction;
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
    T xi=(T)6*sqrt((T)TV::m)*example.grid.One_Over_DX().Min();
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
//#####################################################################
// Function Approximate_Exponential
//#####################################################################
template<class TV> MATRIX<typename TV::SCALAR,TV::m> MPM_DRIVER<TV>::
Approximate_Exponential(const MATRIX<T,TV::m>& A)
{
    MATRIX<T,TV::m> E=A+1;
    if(E.Determinant()>0) return E;
    MATRIX<T,TV::m> G=Approximate_Exponential((T).5*A);
    return G*G;
}
//#####################################################################
namespace PhysBAM{
template class MPM_DRIVER<VECTOR<float,2> >;
template class MPM_DRIVER<VECTOR<float,3> >;
template class MPM_DRIVER<VECTOR<double,2> >;
template class MPM_DRIVER<VECTOR<double,3> >;
}
