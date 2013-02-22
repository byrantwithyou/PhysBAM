//#####################################################################
// Copyright 2013, Chenfanfu Jiang
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Krylov_Solvers/CONJUGATE_GRADIENT.h>
#include <PhysBAM_Tools/Krylov_Solvers/CONJUGATE_RESIDUAL.h>
#include <PhysBAM_Tools/Matrices/DIAGONAL_MATRIX.h>
#include <PhysBAM_Tools/Matrices/MATRIX.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX.h>
#include <PhysBAM_Tools/Read_Write/OCTAVE_OUTPUT.h>
#include <PhysBAM_Tools/Utilities/DEBUG_CAST.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include "MPM_SIMULATION.h"
#include "MPM_SYSTEM.h"
#include "MPM_VECTOR.h"
#include "TIMING.h"
#include <omp.h>
namespace PhysBAM{
using ::std::exp;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> MPM_SIMULATION<TV>::
MPM_SIMULATION()
    :dump_matrix(false),test_system(false),min_mass(1e-8),min_rho((T)0),use_plasticity_yield(true),use_gravity(true),FLIP_alpha((T)0.95),friction_coefficient((T)0.6)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> MPM_SIMULATION<TV>::
~MPM_SIMULATION()
{}
//#####################################################################
// Function Initialize
//#####################################################################
template<class TV> void MPM_SIMULATION<TV>::
Initialize()
{
    LOG::cout<<"Allocating Memories for Simulation..."<<std::endl;
    TIMING_START;
    gravity_constant=TV();gravity_constant(1)=-(T)0.8;
    mu.Resize(particles.number);
    lambda.Resize(particles.number);
    Je.Resize(particles.number);
    Jp.Resize(particles.number);
    Ue.Resize(particles.number);
    Ve.Resize(particles.number);
    Re.Resize(particles.number);
    Se.Resize(particles.number);
    SIGMAe.Resize(particles.number);
    influence_corner.Resize(particles.number);
    weight.Resize(particles.number);
    grad_weight.Resize(particles.number);
    for(int p=0;p<particles.number;p++){
        weight(p).Resize(RANGE<TV_INT>(TV_INT(),TV_INT()+IN));
        grad_weight(p).Resize(RANGE<TV_INT>(TV_INT(),TV_INT()+IN));}
    node_mass.Resize(RANGE<TV_INT>(TV_INT(),grid.counts));
    node_V.Resize(RANGE<TV_INT>(TV_INT(),grid.counts));
    node_V_star.Resize(RANGE<TV_INT>(TV_INT(),grid.counts));
    node_V_old.Resize(RANGE<TV_INT>(TV_INT(),grid.counts));
    node_force.Resize(RANGE<TV_INT>(TV_INT(),grid.counts));
    node_external_force.Resize(RANGE<TV_INT>(TV_INT(),grid.counts));
    frame=0;
    min_mass=particles.mass.Min()*(T)1e-5;
    mu0=ym0/((T)2*((T)1+pr0));
    lambda0=ym0*pr0/(((T)1+pr0)*((T)1-2*pr0));
    if(PROFILING) TIMING_END("");
}
//#####################################################################
// Function Advance_One_Time_Step_Forward_Euler
//#####################################################################
template<class TV> void MPM_SIMULATION<TV>::
Advance_One_Time_Step_Forward_Euler()
{
    Build_Weights_And_Grad_Weights();
    Build_Helper_Structures_For_Constitutive_Model();
    Rasterize_Particle_Data_To_The_Grid();
    if(frame==0) Compute_Particle_Volumes_And_Densities();
    Compute_Grid_Forces();
    if(use_gravity) Apply_Gravity_To_Grid_Forces();
    Update_Velocities_On_Grid();
    Grid_Based_Body_Collisions();
    Solve_The_Linear_System_Explicit();
    T max_node_v=Get_Maximum_Node_Velocity();
    LOG::cout<<"Maximum node velocity: "<<max_node_v<<" = "<<max_node_v/(grid.dX(0)/dt)<<" h/dt"<<std::endl;
    Update_Deformation_Gradient();
    Update_Particle_Velocities();
    T max_particle_v=Get_Maximum_Particle_Velocity();
    LOG::cout<<"Maximum particle velocity: "<<max_particle_v<<" = "<<max_particle_v/(grid.dX.Min()/dt)<<" h/dt"<<std::endl;
    Particle_Based_Body_Collisions();
    Update_Particle_Positions();
    Update_Dirichlet_Box_Positions();
    Update_Colliding_Object_Positions();
    frame++;
}
//#####################################################################
// Function Advance_One_Time_Step_Backward_Euler
//#####################################################################
template<class TV> void MPM_SIMULATION<TV>::
Advance_One_Time_Step_Backward_Euler()
{
    Build_Weights_And_Grad_Weights();
    Build_Helper_Structures_For_Constitutive_Model();
    Rasterize_Particle_Data_To_The_Grid();
    if(frame==0) Compute_Particle_Volumes_And_Densities();
    Compute_Grid_Forces();
    if(use_gravity) Apply_Gravity_To_Grid_Forces();
    Update_Velocities_On_Grid();
    Grid_Based_Body_Collisions();
    Solve_The_Linear_System();
    T max_node_v=Get_Maximum_Node_Velocity();
    LOG::cout<<"Maximum node velocity: "<<max_node_v<<" = "<<max_node_v/(grid.dX(0)/dt)<<" h/dt"<<std::endl;
    Update_Deformation_Gradient();
    Update_Particle_Velocities();
    T max_particle_v=Get_Maximum_Particle_Velocity();
    LOG::cout<<"Maximum particle velocity: "<<max_particle_v<<" = "<<max_particle_v/(grid.dX.Min()/dt)<<" h/dt"<<std::endl;
    Particle_Based_Body_Collisions();
    Update_Particle_Positions();
    Update_Dirichlet_Box_Positions();
    Update_Colliding_Object_Positions();
    frame++;
}
//#####################################################################
// Function Build_Weights_And_Grad_Weights
//#####################################################################
template<class TV> void MPM_SIMULATION<TV>::
Build_Weights_And_Grad_Weights()
{
    TIMING_START;
#pragma omp parallel for
    for(int p=0;p<particles.number;p++){
        grid_basis_function.Build_Weights_And_Grad_Weights_Exact(particles.X(p),grid,influence_corner(p),weight(p),grad_weight(p));}
    if(PROFILING) TIMING_END("Build_Weights_And_Grad_Weights");
}
//#####################################################################
// Function Build_Helper_Structures_For_Constitutive_Model
//#####################################################################
template<class TV> void MPM_SIMULATION<TV>::
Build_Helper_Structures_For_Constitutive_Model()
{
    TIMING_START;
#pragma omp parallel for
    for(int p=0;p<particles.number;p++){
        constitutive_model.Compute_Helper_Quantities_Using_F(particles.Fe(p),particles.Fp(p),Je(p),Jp(p),Ue(p),SIGMAe(p),Ve(p),Re(p),Se(p));
        T lame_scale=exp(xi*((T)1-Jp(p)));
        mu(p)=mu0*lame_scale;
        lambda(p)=lambda0*lame_scale;}
    if(PROFILING) TIMING_END("Build_Helper_Structures_For_Constitutive_Model");
}
//#####################################################################
// Function Rasterize_Particle_Data_To_The_Grid
//#####################################################################
template<class TV> void MPM_SIMULATION<TV>::
Rasterize_Particle_Data_To_The_Grid()
{
    TIMING_START;
    node_mass.Fill(T(0));
    for(int p=0;p<particles.number;p++){
        for(RANGE_ITERATOR<TV::m> it(RANGE<TV_INT>(TV_INT(),TV_INT()+IN));it.Valid();it.Next()){
            node_mass(influence_corner(p)+it.index)+=particles.mass(p)*weight(p)(it.index);}}
    //DEBUG check mass conservation
    // LOG::cout<<"[DEBUG] mass difference grid and particles: "<<node_mass.array.Sum()<<"-"<<particles.mass.Sum()<<"="<<node_mass.array.Sum()-particles.mass.Sum()<<std::endl;
    node_V.Fill(TV());
    for(int p=0;p<particles.number;p++){
        for(RANGE_ITERATOR<TV::m> it(RANGE<TV_INT>(TV_INT(),TV_INT()+IN));it.Valid();it.Next()){
            TV_INT ind=influence_corner(p)+it.index;
            if(node_mass(ind)>min_mass) node_V(ind)+=particles.V(p)*particles.mass(p)*weight(p)(it.index)/node_mass(ind);}}
    if(PROFILING) TIMING_END("Rasterize_Particle_Data_To_The_Grid");
}
//#####################################################################
// Function Compute_Particle_Volumes_And_Densities
//#####################################################################
template<class TV> void MPM_SIMULATION<TV>::
Compute_Particle_Volumes_And_Densities()
{
    TIMING_START;
    T one_over_cell_volume=grid.one_over_dX.Product();
    ARRAY<T> particles_density(particles.number);
    particles_density.Fill(T(0));
    particles.volume.Fill(T(0));
    for(int p=0;p<particles.number;p++){
        for(RANGE_ITERATOR<TV::m> it(RANGE<TV_INT>(TV_INT(),TV_INT()+IN));it.Valid();it.Next()){
            TV_INT ind=influence_corner(p)+it.index;
            particles_density(p)+=node_mass(ind)*weight(p)(it.index)*one_over_cell_volume;}
        if(particles_density(p)>min_rho) particles.volume(p)=particles.mass(p)/particles_density(p);}
    if(PROFILING) TIMING_END("Compute_Particle_Volumes_And_Densities");
}
//#####################################################################
// Function Compute_Grid_Forces
//#####################################################################
template<class TV> void MPM_SIMULATION<TV>::
Compute_Grid_Forces()
{
    TIMING_START;
    node_force.Fill(TV());
    for(int p=0;p<particles.number;p++){
        MATRIX<T,TV::m> B=particles.volume(p)*constitutive_model.Compute_dPsi_dFe(mu(p),lambda(p),particles.Fe(p),Re(p),Je(p)).Times_Transpose(particles.Fe(p));
        for(RANGE_ITERATOR<TV::m> it(RANGE<TV_INT>(TV_INT(),TV_INT()+IN));it.Valid();it.Next()){
            TV_INT ind=influence_corner(p)+it.index;
            node_force(ind)-=B*grad_weight(p)(it.index);}}
    if(PROFILING) TIMING_END("Compute_Grid_Forces");
}
//#####################################################################
// Function Apply_Gravity_To_Grid_Forces
//#####################################################################
template<class TV> void MPM_SIMULATION<TV>::
Apply_Gravity_To_Grid_Forces()
{
    TIMING_START;
    node_external_force.Fill(TV());
    for(int p=0;p<particles.number;p++){
        for(RANGE_ITERATOR<TV::m> it(RANGE<TV_INT>(TV_INT(),TV_INT()+IN));it.Valid();it.Next()){
            TV_INT ind=influence_corner(p)+it.index;
            node_external_force(ind)=node_mass(ind)*gravity_constant;}}
    node_force+=node_external_force;
    if(PROFILING) TIMING_END("Apply_Gravity_To_Grid_Forces");
}
//#####################################################################
// Function Update_Velocities_On_Grid
//#####################################################################
template<class TV> void MPM_SIMULATION<TV>::
Update_Velocities_On_Grid()
{
    TIMING_START;
    for(RANGE_ITERATOR<TV::m> it(RANGE<TV_INT>(TV_INT(),grid.counts));it.Valid();it.Next()){
        node_V_star(it.index)=node_V(it.index);
        if(node_mass(it.index)>min_mass){
            node_V_star(it.index)+=dt/node_mass(it.index)*node_force(it.index);
            for(int b=0;b<dirichlet_box.m;b++)
                if(dirichlet_box(b).Lazy_Inside(grid.Node(it.index)))
                    node_V_star(it.index)=dirichlet_velocity(b);}}
    if(PROFILING) TIMING_END("Update_Velocities_On_Grid");
}
//#####################################################################
// Function Grid_Based_Body_Collisions
//#####################################################################
template<class TV> void MPM_SIMULATION<TV>::
Grid_Based_Body_Collisions()
{
    TIMING_START;
    static T eps=1e-8;
    for(RANGE_ITERATOR<TV::m> it(RANGE<TV_INT>(TV_INT(),grid.counts));it.Valid();it.Next()){
        if(node_mass(it.index)<min_mass) continue;
        const TV& x=grid.Node(it.index);
        if(x(1)<=ground_level && node_V_star(it.index)(1)<=(T)0){
            TV vt(node_V_star(it.index));vt(1)=(T)0;
            T vn=node_V_star(it.index)(1);
            T vt_mag=vt.Magnitude();
            if(vt_mag>eps) node_V_star(it.index)=vt+friction_coefficient*vn*vt/vt_mag;
            else node_V_star(it.index)=vt;}
        for(int b=0;b<rigid_ball.m;b++){
            if(rigid_ball(b).Lazy_Inside(x)){
                TV n=rigid_ball(b).Normal(x);
                TV v_rel=node_V_star(it.index)-rigid_ball_velocity(b);
                T vn=TV::Dot_Product(v_rel,n);
                if(vn<(T)0){
                    TV vt=v_rel-n*vn;
                    T vt_mag=vt.Magnitude();
                    if(vt_mag>eps) node_V_star(it.index)=vt+friction_coefficient*vn*vt/vt_mag+rigid_ball_velocity(b);
                    else node_V_star(it.index)=vt+rigid_ball_velocity(b);}}}}
    if(PROFILING) TIMING_END("Grid_Based_Body_Collisions");
}
//#####################################################################
// Function Solve_The_Linear_System_Explicit
//#####################################################################
template<class TV> void MPM_SIMULATION<TV>::
Solve_The_Linear_System_Explicit()
{
    TIMING_START;
    node_V_old=node_V;node_V=node_V_star;
    if(PROFILING) TIMING_END("Solve_The_Linear_System_Explicit");
}
//#####################################################################
// Function Solve_The_Linear_System
//#####################################################################
template<class TV> void MPM_SIMULATION<TV>::
Solve_The_Linear_System()
{
    TIMING_START;
    static int solve_id=-1;solve_id++;
    MPM_SYSTEM<TV> system(debug_cast<MPM_SIMULATION<TV>&>(*this));
    MPM_VECTOR<TV> rhs,x;
    ARRAY<KRYLOV_VECTOR_BASE<T>*> vectors;
    rhs.v.Resize(RANGE<TV_INT>(TV_INT(),grid.counts));
    x.v.Resize(RANGE<TV_INT>(TV_INT(),grid.counts));
    KRYLOV_SOLVER<T>::Ensure_Size(vectors,x,3);
    rhs.v=node_V_star;
    if(test_system) system.Test_System(*vectors(0),*vectors(1),*vectors(2));
    CONJUGATE_GRADIENT<T> cg;
    CONJUGATE_RESIDUAL<T> cr;
    KRYLOV_SOLVER<T>* solver=&cg;
    solver->print_residuals=true;

    if(dump_matrix){
        LOG::cout<<"solve id "<<solve_id<<std::endl;
        OCTAVE_OUTPUT<T>(STRING_UTILITIES::string_sprintf("M-%i.txt",solve_id).c_str()).Write("M",system,*vectors(0),*vectors(1));
        OCTAVE_OUTPUT<T>(STRING_UTILITIES::string_sprintf("m-%i.txt",solve_id).c_str()).Write("m",node_mass.array);
        OCTAVE_OUTPUT<T>(STRING_UTILITIES::string_sprintf("b-%i.txt",solve_id).c_str()).Write("b",rhs);}

    solver->Solve(system,x,rhs,vectors,(T)1e-7,0,1000);

    if(dump_matrix){
        LOG::cout<<"solve id "<<solve_id<<std::endl;
        OCTAVE_OUTPUT<T>(STRING_UTILITIES::string_sprintf("x-%i.txt",solve_id).c_str()).Write("x",x);}

    node_V_old=node_V;
    node_V=x.v;
    if(PROFILING) TIMING_END("Solve_The_Linear_System");
}
//#####################################################################
// Function Update_Deformation_Gradient
//#####################################################################
template<class TV> void MPM_SIMULATION<TV>::
Update_Deformation_Gradient()
{
    TIMING_START;
    if(!use_plasticity_yield){
#pragma omp parallel for
        for(int p=0;p<particles.number;p++){
            MATRIX<T,TV::m> grad_vp;
            for(RANGE_ITERATOR<TV::m> it(RANGE<TV_INT>(TV_INT(),TV_INT()+IN));it.Valid();it.Next()){
                TV_INT ind=influence_corner(p)+it.index;
                grad_vp+=MATRIX<T,TV::m>::Outer_Product(node_V(ind),grad_weight(p)(it.index));}
            particles.Fe(p)=particles.Fe(p)+dt*grad_vp*particles.Fe(p);}}
    else{
#pragma omp parallel for
        for(int p=0;p<particles.number;p++){
            MATRIX<T,TV::m> grad_vp;
            for(RANGE_ITERATOR<TV::m> it(RANGE<TV_INT>(TV_INT(),TV_INT()+IN));it.Valid();it.Next()){
                TV_INT ind=influence_corner(p)+it.index;
                grad_vp+=MATRIX<T,TV::m>::Outer_Product(node_V(ind),grad_weight(p)(it.index));}
            MATRIX<T,TV::m> Fe_hat=particles.Fe(p)+dt*grad_vp*particles.Fe(p);
            MATRIX<T,TV::m> Fp_hat=particles.Fp(p);
            MATRIX<T,TV::m> F=Fe_hat*Fp_hat;
            MATRIX<T,TV::m> U_hat,V_hat;
            DIAGONAL_MATRIX<T,TV::m> SIGMA_hat;
            Fe_hat.Fast_Singular_Value_Decomposition(U_hat,SIGMA_hat,V_hat);
            SIGMA_hat=SIGMA_hat.Clamp_Min(yield_min);
            SIGMA_hat=SIGMA_hat.Clamp_Max(yield_max);
            particles.Fe(p)=U_hat*SIGMA_hat.Times_Transpose(V_hat);
            particles.Fp(p)=V_hat*(SIGMA_hat.Inverse()).Times_Transpose(U_hat)*F;
            // if(use_plasticity_clamp){
            //     MATRIX<T,TV::m> Uphat,Vphat;
            //     DIAGONAL_MATRIX<T,TV::m> SIGMAphat;
            //     particles.Fp(p).Fast_Singular_Value_Decomposition(Uphat,SIGMAphat,Vphat);
            //     if(SIGMAphat.Min()<clamp_min || SIGMAphat.Max()>clamp_max){
            //         particles.Fp(p)=Fp_hat;
            //         particles.Fe(p)=Fe_hat;}}
            if(use_plasticity_clamp){
                MATRIX<T,TV::m> Uphat,Vphat;
                DIAGONAL_MATRIX<T,TV::m> SIGMAphat;
                particles.Fp(p).Fast_Singular_Value_Decomposition(Uphat,SIGMAphat,Vphat);
                SIGMAphat.Clamp_Min(clamp_min);
                SIGMAphat.Clamp_Max(clamp_max);
                particles.Fp(p)=Uphat*SIGMAphat.Times_Transpose(Vphat);
                particles.Fe(p)=F*(particles.Fp(p).Inverse());}
        }}
    if(PROFILING) TIMING_END("Update_Deformation_Gradient");
}
//#####################################################################
// Function Update_Particle_Velocities
//#####################################################################
template<class TV> void MPM_SIMULATION<TV>::
Update_Particle_Velocities()
{
    TIMING_START;
#pragma omp parallel for
    for(int p=0;p<particles.number;p++){
        TV V_PIC;
        TV V_FLIP=particles.V(p);
        for(RANGE_ITERATOR<TV::m> it(RANGE<TV_INT>(TV_INT(),TV_INT()+IN));it.Valid();it.Next()){
            TV_INT ind=influence_corner(p)+it.index;
            T w=weight(p)(it.index);
            V_PIC+=node_V(ind)*w;
            V_FLIP+=(node_V(ind)-node_V_old(ind))*w;}
        particles.V(p)=((T)1-FLIP_alpha)*V_PIC+FLIP_alpha*V_FLIP;}
    if(PROFILING) TIMING_END("Update_Particle_Velocities");
}
//#####################################################################
// Function Particle_Based_Body_Collisions
//#####################################################################
template<class TV> void MPM_SIMULATION<TV>::
Particle_Based_Body_Collisions()
{
    TIMING_START;
    static T eps=1e-8;
#pragma omp parallel for
    for(int p=0;p<particles.number;p++){
        TV& x=particles.X(p);
        if(x(1)<=ground_level && particles.V(p)(1)<=(T)0){
            TV vt(particles.V(p));vt(1)=(T)0;
            T vn=particles.V(p)(1);
            T vt_mag=vt.Magnitude();
            if(vt_mag>eps) particles.V(p)=vt+friction_coefficient*vn*vt/vt_mag;
            else particles.V(p)=vt;}
        for(int b=0;b<rigid_ball.m;b++){
            if(rigid_ball(b).Lazy_Inside(x)){
                TV n=rigid_ball(b).Normal(x);
                TV v_rel=particles.V(p)-rigid_ball_velocity(b);
                T vn=TV::Dot_Product(v_rel,n);
                if(vn<(T)0){
                    TV vt=v_rel-n*vn;
                    T vt_mag=vt.Magnitude();
                    if(vt_mag>eps) particles.V(p)=vt+friction_coefficient*vn*vt/vt_mag+rigid_ball_velocity(b);
                    else particles.V(p)=vt+rigid_ball_velocity(b);}}}}
    if(PROFILING) TIMING_END("Particle_Based_Body_Collisions");
}
//#####################################################################
// Function Update_Particle_Positions
//#####################################################################
template<class TV> void MPM_SIMULATION<TV>::
Update_Particle_Positions()
{
    TIMING_START;
#pragma omp parallel for
    for(int p=0;p<particles.number;p++)
        particles.X(p)+=dt*particles.V(p);
    if(PROFILING) TIMING_END("Update_Particle_Positions");
}
//#####################################################################
// Function Update_Dirichlet_Box_Positions
//#####################################################################
template<class TV> void MPM_SIMULATION<TV>::
Update_Dirichlet_Box_Positions()
{
    TIMING_START;
    for(int b=0;b<dirichlet_box.m;b++)
        dirichlet_box(b)+=dt*dirichlet_velocity(b);
    if(PROFILING) TIMING_END("Update_Particle_Positions");
}
//#####################################################################
// Function Update_Colliding_Object_Positions
//#####################################################################
template<class TV> void MPM_SIMULATION<TV>::
Update_Colliding_Object_Positions()
{
    TIMING_START;
    for(int b=0;b<rigid_ball.m;b++)
        rigid_ball(b).center+=dt*rigid_ball_velocity(b);
    if(PROFILING) TIMING_END("Update_Particle_Positions");
}
//#####################################################################
// Function Get_Maximum_Node_Velocity
//#####################################################################
template<class TV> typename TV::SCALAR MPM_SIMULATION<TV>::
Get_Maximum_Node_Velocity() const
{
    T max_node_v=0;
    for(int i=0;i<node_V.array.m;i++){
        T mag=node_V.array(i).Magnitude();
        if(mag>max_node_v) max_node_v=mag;}
    return max_node_v;
}
//#####################################################################
// Function Get_Maximum_Particle_Velocity
//#####################################################################
template<class TV> typename TV::SCALAR MPM_SIMULATION<TV>::
Get_Maximum_Particle_Velocity() const
{
    T max_particle_v=0;
    for(int i=0;i<particles.V.m;i++){
        T mag=particles.V(i).Magnitude();
        if(mag>max_particle_v) max_particle_v=mag;}
    return max_particle_v;
}
//#####################################################################
template class MPM_SIMULATION<VECTOR<float,2> >;
template class MPM_SIMULATION<VECTOR<float,3> >;
template class MPM_SIMULATION<VECTOR<double,2> >;
template class MPM_SIMULATION<VECTOR<double,3> >;
}
