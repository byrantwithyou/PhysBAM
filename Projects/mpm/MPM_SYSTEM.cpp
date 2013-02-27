//#####################################################################
// Copyright 2013, Chenfanfu Jiang
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <PhysBAM_Tools/Utilities/DEBUG_CAST.h>
#include <iomanip>
#include "MPM_SIMULATION.h"
#include "MPM_SYSTEM.h"
#include "MPM_VECTOR.h"
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> MPM_SYSTEM<TV>::
MPM_SYSTEM(const MPM_SIMULATION<TV>& sim)
    :BASE(false,true),sim(sim),beta(1)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> MPM_SYSTEM<TV>::
~MPM_SYSTEM()
{
}
//#####################################################################
// Function Force
//#####################################################################
template<class TV> void MPM_SYSTEM<TV>::
Force(ARRAY<TV,TV_INT>& f) const
{
    f.Fill(TV());
    for(int p=0;p<sim.particles.number;p++){
        MATRIX<T,TV::m> B=sim.particles.volume(p)*sim.constitutive_model.Compute_dPsi_dFe(sim.mu(p),sim.lambda(p),sim.particles.Fe(p),sim.Re(p),sim.Je(p))*(sim.particles.Fe(p).Transposed());
        for(RANGE_ITERATOR<TV::m> it(RANGE<TV_INT>(TV_INT(),TV_INT()+sim.IN));it.Valid();it.Next())
            f(sim.influence_corner(p)+it.index)-=B*sim.grad_weight(p)(it.index);}
}
//#####################################################################
// Function Apply_Force_Derivatives
//#####################################################################
template<class TV> void MPM_SYSTEM<TV>::
Apply_Force_Derivatives(const ARRAY<TV,TV_INT>& du,ARRAY<TV,TV_INT>& df) const
{ 
    df.Resize(RANGE<TV_INT>(TV_INT(),sim.grid.counts));
    df.Fill(TV());
    for(int p=0;p<sim.particles.number;p++){
        MATRIX<T,TV::m> Cp;
        for(RANGE_ITERATOR<TV::m> it(RANGE<TV_INT>(TV_INT(),TV_INT()+sim.IN));it.Valid();it.Next()){
            bool contribute=true;
            // for(int b=0;b<sim.dirichlet_box.m;b++) if(sim.dirichlet_box(b).Lazy_Inside(sim.grid.Node(sim.influence_corner(p)+it.index))) contribute=false;
            // for(int b=0;b<sim.rigid_ball.m;b++) if(sim.rigid_ball(b).Lazy_Inside(sim.grid.Node(sim.influence_corner(p)+it.index))) contribute=false;
            // for(int b=0;b<sim.rigid_box.m;b++) if(sim.rigid_box(b).Lazy_Inside(sim.grid.Node(sim.influence_corner(p)+it.index))) contribute=false;
            // if(sim.grid.Node(sim.influence_corner(p)+it.index)(1)<sim.ground_level) contribute=false;
            if(contribute) Cp+=MATRIX<T,TV::m>::Outer_Product(du(sim.influence_corner(p)+it.index),sim.grad_weight(p)(it.index));}
        MATRIX<T,TV::m> Ep=Cp*sim.particles.Fe(p);
        MATRIX<T,TV::m> Ap=sim.constitutive_model.Compute_d2Psi_dFe_dFe_Action_dF(sim.mu(p),sim.lambda(p),sim.particles.Fe(p),sim.Je(p),sim.Re(p),sim.Se(p),Ep);
        MATRIX<T,TV::m> Gp=sim.particles.volume(p)*(Ap.Times_Transpose(sim.particles.Fe(p)));
        for(RANGE_ITERATOR<TV::m> it(RANGE<TV_INT>(TV_INT(),TV_INT()+sim.IN));it.Valid();it.Next()){
            bool contribute=true;
            // for(int b=0;b<sim.dirichlet_box.m;b++) if(sim.dirichlet_box(b).Lazy_Inside(sim.grid.Node(sim.influence_corner(p)+it.index))) contribute=false;
            // for(int b=0;b<sim.rigid_ball.m;b++) if(sim.rigid_ball(b).Lazy_Inside(sim.grid.Node(sim.influence_corner(p)+it.index))) contribute=false;
            // for(int b=0;b<sim.rigid_box.m;b++) if(sim.rigid_box(b).Lazy_Inside(sim.grid.Node(sim.influence_corner(p)+it.index))) contribute=false;
            // if(sim.grid.Node(sim.influence_corner(p)+it.index)(1)<sim.ground_level) contribute=false;
            if(contribute) df(sim.influence_corner(p)+it.index)-=Gp*sim.grad_weight(p)(it.index);}}
}
//#####################################################################
// Function Multiply
//#####################################################################
template<class TV> void MPM_SYSTEM<TV>::
Multiply(const KRYLOV_VECTOR_BASE<T>& x,KRYLOV_VECTOR_BASE<T>& result) const
{
    ARRAY<TV,TV_INT>& rr=debug_cast<MPM_VECTOR<TV>&>(result).v;
    const ARRAY<TV,TV_INT>& xx=debug_cast<const MPM_VECTOR<TV>&>(x).v;
    Apply_Force_Derivatives(xx,rr);
    T beta_dt2=beta*sqr(sim.dt);
    for(RANGE_ITERATOR<TV::m> it(RANGE<TV_INT>(TV_INT(),TV_INT()+sim.grid.counts));it.Valid();it.Next()){
        if(sim.node_mass(it.index)>sim.min_mass)
            rr(it.index)=xx(it.index)+beta_dt2/sim.node_mass(it.index)*rr(it.index);
        else rr(it.index)=TV();}
}
//#####################################################################
// Function Project
//#####################################################################
template<class TV> void MPM_SYSTEM<TV>::
Project(KRYLOV_VECTOR_BASE<T>& x) const
{
    ARRAY<TV,TV_INT>& xx=debug_cast<MPM_VECTOR<TV>&>(x).v;
    for(RANGE_ITERATOR<TV::m> it(RANGE<TV_INT>(TV_INT(),TV_INT()+sim.grid.counts));it.Valid();it.Next())
        for(int b=0;b<sim.dirichlet_box.m;b++)
            if(sim.dirichlet_box(b).Lazy_Inside(sim.grid.Node(it.index)))
                xx(it.index)=TV();
}
//#####################################################################
// Function Inner_Product
//#####################################################################
template<class TV> double MPM_SYSTEM<TV>::
Inner_Product(const KRYLOV_VECTOR_BASE<T>& x,const KRYLOV_VECTOR_BASE<T>& y) const
{
    const ARRAY<TV,TV_INT>& xx=debug_cast<const MPM_VECTOR<TV>&>(x).v;
    const ARRAY<TV,TV_INT>& yy=debug_cast<const MPM_VECTOR<TV>&>(y).v;
    double r=0;
    for(int i=0;i<xx.array.m;i++) r+=sim.node_mass.array(i)*xx.array(i).Dot(yy.array(i));
    return r;
}
//#####################################################################
// Function Convergence_Norm
//#####################################################################
template<class TV> typename TV::SCALAR MPM_SYSTEM<TV>::
Convergence_Norm(const KRYLOV_VECTOR_BASE<T>& x) const
{
    const ARRAY<TV,TV_INT>& xx=debug_cast<const MPM_VECTOR<TV>&>(x).v;
    return xx.array.Maximum_Magnitude();
}
//#####################################################################
// Function Set_Boundary_Conditions
//#####################################################################
template<class TV> void MPM_SYSTEM<TV>::
Set_Boundary_Conditions(KRYLOV_VECTOR_BASE<T>& x) const
{
}
//#####################################################################
// Function Project_Nullspace
//#####################################################################
template<class TV> void MPM_SYSTEM<TV>::
Project_Nullspace(KRYLOV_VECTOR_BASE<T>& x) const
{
}
//#####################################################################
// Function Apply_Preconditioner
//#####################################################################
template<class TV> void MPM_SYSTEM<TV>::
Apply_Preconditioner(const KRYLOV_VECTOR_BASE<T>& r,KRYLOV_VECTOR_BASE<T>& z) const
{
}
namespace PhysBAM{
template class MPM_SYSTEM<VECTOR<float,2> >;
template class MPM_SYSTEM<VECTOR<float,3> >;
template class MPM_SYSTEM<VECTOR<double,2> >;
template class MPM_SYSTEM<VECTOR<double,3> >;
}
