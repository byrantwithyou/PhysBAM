//#####################################################################
// Copyright 2013, Chenfanfu Jiang
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Krylov_Solvers/CONJUGATE_GRADIENT.h>
#include <PhysBAM_Tools/Krylov_Solvers/CONJUGATE_RESIDUAL.h>
#include <PhysBAM_Tools/Math_Tools/cbrt.h>
#include <PhysBAM_Tools/Math_Tools/clamp.h>
#include <PhysBAM_Tools/Math_Tools/Inverse.h>
#include <PhysBAM_Tools/Math_Tools/pow.h>
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
    :dump_matrix(false),test_system(false),min_mass(1e-8),min_rho((T)0),use_gravity(true),FLIP_alpha((T)0.95),friction_coefficient((T)0.8)
{
    constitutive_model.dev_part_only=true;
}
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
    gravity_constant=TV();gravity_constant(1)=-(T)9.8;
    Je.Resize(particles.number);
    Re.Resize(particles.number);
    Se.Resize(particles.number);
    influence_corner.Resize(particles.number);
    weight.Resize(particles.number);
    grad_weight.Resize(particles.number);
    for(int p=0;p<particles.number;p++){
        weight(p).Resize(RANGE<TV_INT>(TV_INT(),TV_INT()+IN));
        grad_weight(p).Resize(RANGE<TV_INT>(TV_INT(),TV_INT()+IN));}
    node_mass.Resize(RANGE<TV_INT>(TV_INT(),grid.counts));
    node_volume.Resize(RANGE<TV_INT>(TV_INT(),grid.counts));
    node_V.Resize(RANGE<TV_INT>(TV_INT(),grid.counts));
    node_V_star.Resize(RANGE<TV_INT>(TV_INT(),grid.counts));
    node_V_old.Resize(RANGE<TV_INT>(TV_INT(),grid.counts));
    node_force.Resize(RANGE<TV_INT>(TV_INT(),grid.counts));
    frame=0;
    min_mass=particles.mass.Min()*(T)1e-5;
    min_volume=(T)1e-7;
    min_density=(T)1e-7;
}
//#####################################################################
// Function Resize_For_New_Particle_Data
//#####################################################################
template<class TV> void MPM_SIMULATION<TV>::
Resize_For_New_Particle_Data()
{
    Je.Resize(particles.number);
    Re.Resize(particles.number);
    Se.Resize(particles.number);
    influence_corner.Resize(particles.number);
    weight.Resize(particles.number);
    grad_weight.Resize(particles.number);
    for(int p=0;p<particles.number;p++){
        weight(p).Resize(RANGE<TV_INT>(TV_INT(),TV_INT()+IN));
        grad_weight(p).Resize(RANGE<TV_INT>(TV_INT(),TV_INT()+IN));}
    min_mass=particles.mass.Min()*(T)1e-5;
    min_volume=(T)1e-7;
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
        constitutive_model.Compute_Helper_Quantities_Using_F(particles.Fe(p),Je(p),Re(p),Se(p));
        T lame_scale=exp(xi*(1-particles.Fp(p).Determinant()));
        particles.mu(p)=particles.mu0(p)*lame_scale;
        particles.lambda(p)=particles.lambda0(p)*lame_scale;}
}
//#####################################################################
// Function Build_Pressure_And_One_Over_Lambda_J
//#####################################################################
template<class TV> void MPM_SIMULATION<TV>::
Build_Pressure_And_One_Over_Lambda_J()
{
#pragma omp parallel for schedule(static)
    for(int p=0;p<particles.number;p++){
        if(particles.compress(p)){
            T J_lazy=Je(p);
            // particles.pressure(p)=particles.lambda(p)*(J_lazy-1);
            particles.one_over_lambda_J(p)=(T)1.0/(particles.lambda(p)*J_lazy);}
            // LOG::cout<<particles.one_over_lambda_J(p)<<std::endl;}
        
        else{
            // particles.pressure(p)=(T)0;
            particles.one_over_lambda_J(p)=(T)0;}}
 }
    
//#####################################################################
// Function Rasterize_Particle_Data_To_The_Grid
//#####################################################################
template<class TV> void MPM_SIMULATION<TV>::
Rasterize_Particle_Data_To_The_Grid()
{
    TIMING_START;
    node_mass.Fill(T(0));
    node_volume.Fill(T(0));
    for(int p=0;p<particles.number;p++){
        for(RANGE_ITERATOR<TV::m> it(RANGE<TV_INT>(TV_INT(),TV_INT()+IN));it.Valid();it.Next()){
            node_mass(influence_corner(p)+it.index)+=particles.mass(p)*weight(p)(it.index);
            node_volume(influence_corner(p)+it.index)+=(T)1.0*weight(p)(it.index);}}
    //DEBUG check mass conservation
    // LOG::cout<<"[DEBUG] mass difference grid and particles: "<<node_mass.array.Sum()<<"-"<<particles.mass.Sum()<<"="<<node_mass.array.Sum()-particles.mass.Sum()<<std::endl;
    node_V.Fill(TV());
    for(int p=0;p<particles.number;p++){
        for(RANGE_ITERATOR<TV::m> it(RANGE<TV_INT>(TV_INT(),TV_INT()+IN));it.Valid();it.Next()){
            TV_INT ind=influence_corner(p)+it.index;
            if(node_mass(ind)>min_mass) node_V(ind)+=particles.V(p)*particles.mass(p)*weight(p)(it.index)/node_mass(ind);}}
            // if(node_mass(ind)!=0) node_V(ind)+=particles.V(p)*particles.mass(p)*weight(p)(it.index)/node_mass(ind);}}
}
//#####################################################################
// Function Compute_Particle_Volumes_And_Densities
//#####################################################################
template<class TV> void MPM_SIMULATION<TV>::
Compute_Particle_Volumes_And_Densities(int starting_index,int count)
{
    TIMING_START;
    T one_over_cell_volume=grid.one_over_dX.Product();
    ARRAY<T> particles_density(count);
#pragma omp parallel for
    for(int p=starting_index;p<starting_index+count;p++){
        particles_density(p-starting_index)=(T)0;
        particles.volume(p)=(T)0;}
#pragma omp parallel for
    for(int p=starting_index;p<starting_index+count;p++){
        for(RANGE_ITERATOR<TV::m> it(RANGE<TV_INT>(TV_INT(),TV_INT()+IN));it.Valid();it.Next()){
            TV_INT ind=influence_corner(p)+it.index;
            particles_density(p-starting_index)+=node_mass(ind)*weight(p)(it.index)*one_over_cell_volume;}
        if(particles_density(p-starting_index)>min_rho) particles.volume(p)=particles.mass(p)/particles_density(p-starting_index);
        else PHYSBAM_FATAL_ERROR();}
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
        MATRIX<T,TV::m> B=particles.volume(p)*constitutive_model.Compute_dPsi_dFe(particles.mu(p),particles.lambda(p),particles.Fe(p),Re(p),Je(p)).Times_Transpose(particles.Fe(p));
        for(RANGE_ITERATOR<TV::m> it(RANGE<TV_INT>(TV_INT(),TV_INT()+IN));it.Valid();it.Next()){
            TV_INT ind=influence_corner(p)+it.index;
            node_force(ind)-=B*grad_weight(p)(it.index);}}
}
//#####################################################################
// Function Apply_Gravity_To_Grid_Forces
//#####################################################################
template<class TV> void MPM_SIMULATION<TV>::
Apply_Gravity_To_Grid_Forces()
{
    TIMING_START;
    for(RANGE_ITERATOR<TV::m> it(RANGE<TV_INT>(TV_INT(),grid.counts));it.Valid();it.Next())
        if(node_mass(it.index)>=min_mass)
            node_force(it.index)+=node_mass(it.index)*gravity_constant;
}
//#####################################################################
// Function Update_Velocities_On_Grid
//#####################################################################
template<class TV> void MPM_SIMULATION<TV>::
Update_Velocities_On_Grid()
{
    TIMING_START;
    nodes_need_projection.Clean_Memory();
    for(RANGE_ITERATOR<TV::m> it(RANGE<TV_INT>(TV_INT(),grid.counts));it.Valid();it.Next()){
        node_V_star(it.index)=node_V(it.index);
        if(node_mass(it.index)>min_mass){
            node_V_star(it.index)+=dt/node_mass(it.index)*node_force(it.index);
            for(int b=0;b<dirichlet_box.m;b++)
                if(dirichlet_box(b).Lazy_Inside(grid.Node(it.index))){
                    node_V_star(it.index)=dirichlet_velocity(b);
                    nodes_need_projection.Append(it.index);
                    break;}}}
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
        if(node_mass(it.index)>min_mass){
            const TV& x=grid.Node(it.index);
            for(int d=0;d<TV::m;d++){
                T left_wall=grid.domain.min_corner(d)+4.2*grid.dX.Min(),right_wall=grid.domain.max_corner(d)-4.2*grid.dX.Min();
                if(x(d)<left_wall && node_V_star(it.index)(d)<(T)0){
                    TV vt(node_V_star(it.index));vt(d)=(T)0;
                    T vn=node_V_star(it.index)(d); // <0
                    T vt_mag=vt.Magnitude();
                    if(vt_mag>eps){
                        T impulse=friction_coefficient*vn/vt_mag;
                        if(impulse<-(T)1) impulse=-(T)1;
                        node_V_star(it.index)=vt*((T)1+impulse);}
                    else node_V_star(it.index)=vt;
                    // node_V_star(it.index)=TV(); // sticky
                }
                if(x(d)>right_wall && node_V_star(it.index)(d)>(T)0){
                    TV vt(node_V_star(it.index));vt(d)=(T)0;
                    T vn=node_V_star(it.index)(d); // >0
                    T vt_mag=vt.Magnitude();
                    if(vt_mag>eps){
                        T impulse=-friction_coefficient*vn/vt_mag;
                        if(impulse<-(T)1) impulse=-(T)1;
                        node_V_star(it.index)=vt*((T)1+impulse);}
                    else node_V_star(it.index)=vt;
                    // node_V_star(it.index)=TV(); // sticky
                }
                nodes_need_projection.Append(it.index);}
            for(int b=0;b<rigid_ball.m;b++){
                if(rigid_ball(b).Lazy_Inside(x)){
                    TV n=rigid_ball(b).Normal(x);
                    TV v_rel=node_V_star(it.index)-rigid_ball_velocity(b);
                    T vn=TV::Dot_Product(v_rel,n);
                    if(vn<(T)0){
                        TV vt=v_rel-n*vn;
                        T vt_mag=vt.Magnitude();
                        if(vt_mag>eps){
                            T impulse=friction_coefficient*vn/vt_mag;
                            if(impulse<-(T)1) impulse=-(T)1;
                            node_V_star(it.index)=vt*((T)1+impulse)+rigid_ball_velocity(b);}
                        else node_V_star(it.index)=vt+rigid_ball_velocity(b);
                        nodes_need_projection.Append(it.index);
                        break;}}}}}
}
//#####################################################################
// Function Solve_The_Linear_System_Explicit
//#####################################################################
template<class TV> void MPM_SIMULATION<TV>::
Solve_The_Linear_System_Explicit()
{
    TIMING_START;
    node_V_old=node_V;node_V=node_V_star;
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
    x.v=rhs.v;
    if(test_system) system.Test_System(*vectors(0),*vectors(1),*vectors(2));
    CONJUGATE_GRADIENT<T> cg;
    CONJUGATE_RESIDUAL<T> cr;
    KRYLOV_SOLVER<T>* solver=&cg;
    solver->print_residuals=false;
    if(dump_matrix){ // load M-1.txt;load m-1.txt;mm=reshape([m m]',rows(M),1);R=diag(mm)*M;max(max(abs(R-R')))
        LOG::cout<<"solve id "<<solve_id<<std::endl;
        OCTAVE_OUTPUT<T>(STRING_UTILITIES::string_sprintf("M-%i.txt",solve_id).c_str()).Write("M",system,*vectors(0),*vectors(1));
        OCTAVE_OUTPUT<T>(STRING_UTILITIES::string_sprintf("m-%i.txt",solve_id).c_str()).Write("m",node_mass.array);
        OCTAVE_OUTPUT<T>(STRING_UTILITIES::string_sprintf("b-%i.txt",solve_id).c_str()).Write("b",rhs);} 
    solver->Solve(system,x,rhs,vectors,(T)1e-12,0,1000);
    if(dump_matrix){
        LOG::cout<<"solve id "<<solve_id<<std::endl;
        OCTAVE_OUTPUT<T>(STRING_UTILITIES::string_sprintf("x-%i.txt",solve_id).c_str()).Write("x",x);}
    node_V_old=node_V;
    node_V=x.v;
    vectors.Delete_Pointers_And_Clean_Memory();
}
//#####################################################################
// Function Update_Deformation_Gradient
//#####################################################################
template<class TV> void MPM_SIMULATION<TV>::
Update_Deformation_Gradient()
{
    TIMING_START;
#pragma omp parallel for
     for(int p=0;p<particles.number;p++){
         if(particles.use_visco_plasticity(p)){
             MATRIX<T,TV::m> grad_vp;
             for(RANGE_ITERATOR<TV::m> it(RANGE<TV_INT>(TV_INT(),TV_INT()+IN));it.Valid();it.Next()){
                 TV_INT ind=influence_corner(p)+it.index;
                 grad_vp+=MATRIX<T,TV::m>::Outer_Product(node_V(ind),grad_weight(p)(it.index));}
             MATRIX<T,TV::m> Fe_hat=particles.Fe(p)+dt*grad_vp*particles.Fe(p);
             MATRIX<T,TV::m> Fp_hat=particles.Fp(p);
             MATRIX<T,TV::m> F=Fe_hat*Fp_hat;
             MATRIX<T,TV::m> U_hat,V_hat;
             DIAGONAL_MATRIX<T,TV::m> Sigma;
             Fe_hat.Fast_Singular_Value_Decomposition(U_hat,Sigma,V_hat);
             T Pnorm=(constitutive_model.Compute_dPsi_dFe(particles.mu(p),particles.lambda(p),particles.Fe(p),Re(p),Je(p))).Frobenius_Norm();
             T gamma=(T)0;
             if(Pnorm>1e-10) gamma=clamp(dt*particles.visco_nu(p)*(Pnorm-particles.visco_tau(p))/Pnorm,(T)0,(T)1);
             T scale=Inverse(cbrt(Sigma.Determinant()));
             DIAGONAL_MATRIX<T,TV::m> Middle;
             for(int d=0;d<TV::m;d++) Middle(d,d)=pow(Sigma(d,d)*scale,gamma);
             particles.Fp(p)=V_hat*Middle.Times_Transpose(V_hat)*Fp_hat;
             particles.Fe(p)=F*particles.Fp(p).Inverse();
             particles.visco_tau(p)+=particles.visco_kappa(p)*gamma*Pnorm;}
         else if(particles.use_plasticity_yield(p)){
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
             SIGMA_hat=SIGMA_hat.Clamp_Min(particles.yield_min(p));
             SIGMA_hat=SIGMA_hat.Clamp_Max(particles.yield_max(p));
             particles.Fe(p)=U_hat*SIGMA_hat.Times_Transpose(V_hat);
             particles.Fp(p)=V_hat*(SIGMA_hat.Inverse()).Times_Transpose(U_hat)*F;
             if(particles.use_plasticity_clamp(p)){
                 MATRIX<T,TV::m> Uphat,Vphat;
                 DIAGONAL_MATRIX<T,TV::m> SIGMAphat;
                 particles.Fp(p).Fast_Singular_Value_Decomposition(Uphat,SIGMAphat,Vphat);
                 SIGMAphat.Clamp_Min(particles.clamp_min(p));
                 SIGMAphat.Clamp_Max(particles.clamp_max(p));
                 particles.Fp(p)=Uphat*SIGMAphat.Times_Transpose(Vphat);
                 particles.Fe(p)=F*(particles.Fp(p).Inverse());}}
         else{
             MATRIX<T,TV::m> grad_vp;grad_vp*=(T)0;
             for(RANGE_ITERATOR<TV::m> it(RANGE<TV_INT>(TV_INT(),TV_INT()+IN));it.Valid();it.Next()){
                 TV_INT ind=influence_corner(p)+it.index;
                 grad_vp+=MATRIX<T,TV::m>::Outer_Product(node_V(ind),grad_weight(p)(it.index));}
             MATRIX<T,TV::m> dFe=dt*grad_vp*particles.Fe(p);
             
             // get rid of numerically small dFe
//             for(int d1=0;d1<TV::m;d1++)
//                 for(int d2=0;d2<TV::m;d2++)
//                     if(abs(dFe(d1,d2))<(T)1e-4)
//                         dFe(d1,d2)=(T)0;
             
             particles.Fe(p)=particles.Fe(p)+dFe;}

     }
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
        TV x=particles.X(p);
        for(int d=0;d<TV::m;d++){
            T left_wall=grid.domain.min_corner(d)+4.2*grid.dX.Min(),right_wall=grid.domain.max_corner(d)-4.2*grid.dX.Min();
            if(x(d)<=left_wall && particles.V(p)(d)>0){
                TV vt(particles.V(p));vt(d)=(T)0;
                particles.V(p)=vt;
            }
            if(x(d)>=right_wall && particles.V(p)(d)<0){
                TV vt(particles.V(p));vt(d)=(T)0;
                particles.V(p)=vt;
            }
        }
        for(int b=0;b<rigid_ball.m;b++){
            if(rigid_ball(b).Lazy_Inside(x)){
                TV n=rigid_ball(b).Normal(x);
                TV v_rel=particles.V(p)-rigid_ball_velocity(b);
                T vn=TV::Dot_Product(v_rel,n);
                if(vn<(T)0){
                    TV vt=v_rel-n*vn;
                    T vt_mag=vt.Magnitude();
                    if(vt_mag>eps){
                        T impulse=friction_coefficient*vn/vt_mag;
                        if(impulse<(T)-1) impulse=(T)-1;
                        particles.V(p)=vt*(1+impulse)+rigid_ball_velocity(b);}
                    else particles.V(p)=vt+rigid_ball_velocity(b);}}}}
}
//#####################################################################
// Function Update_Particle_Positions
//#####################################################################
template<class TV> void MPM_SIMULATION<TV>::
Update_Particle_Positions()
{
    TIMING_START;
#pragma omp parallel for
    for(int p=0;p<particles.number;p++){
        // if(dt*particles.V(p).Max_Mag()>grid.dX.Min()){
        //     LOG::cout<<"breaking CFL"<<std::endl;
        //     PHYSBAM_FATAL_ERROR();}
        particles.X(p)+=dt*particles.V(p);}
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
// Function Get_Total_Momentum_On_Nodes
//#####################################################################
template<class TV> TV MPM_SIMULATION<TV>::
Get_Total_Momentum_On_Nodes() const
{
    TV momen;
    for(int i=0;i<node_V.array.m;i++)
        if(node_mass.array(i)>min_mass)
            momen+=node_mass.array(i)*node_V.array(i);
    return momen;
}
//#####################################################################
// Function Get_Total_Momentum_On_Particles
//#####################################################################
template<class TV> TV MPM_SIMULATION<TV>::
Get_Total_Momentum_On_Particles() const
{
    TV momen;
    for(int i=0;i<particles.V.m;i++)
        momen+=particles.mass(i)*particles.V(i);
    return momen;
}
//#####################################################################
template class MPM_SIMULATION<VECTOR<float,2> >;
template class MPM_SIMULATION<VECTOR<float,3> >;
template class MPM_SIMULATION<VECTOR<double,2> >;
template class MPM_SIMULATION<VECTOR<double,3> >;
}
