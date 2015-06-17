//#####################################################################
// Copyright 2014.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Grids_Uniform/CELL_ITERATOR.h>
#include <Tools/Grids_Uniform_Interpolation/AVERAGING_UNIFORM.h>
#include <Tools/Grids_Uniform_Interpolation/INTERPOLATION_UNIFORM.h>
#include <Tools/Grids_Uniform_Interpolation/QUADRATIC_INTERPOLATION_UNIFORM.h>
#include <Tools/Log/LOG.h>
#include <Tools/Nonlinear_Equations/ITERATIVE_SOLVER.h>
#include <Tools/Nonlinear_Equations/NONLINEAR_FUNCTION.h>
#include <Tools/Vectors/VECTOR.h>
#include <Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <Geometry/Level_Sets/LEVELSET.h>
#include "ADVECTION_UPWIND.h"
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV,class T2,class Z_INTERP,class T_AVERAGING,class U_INTERP,class T_FACE_LOOKUP>
ADVECTION_UPWIND<TV,T2,Z_INTERP,T_AVERAGING,U_INTERP,T_FACE_LOOKUP>::
ADVECTION_UPWIND(const LEVELSET<TV>& levelset,const T_FACE_LOOKUP& face_velocities,
    T max_in,T max_out,std::function<T2(const TV& X,T time)> bc_Z,T time)
    :levelset(levelset),face_velocities(face_velocities),max_in(max_in),max_out(max_out),time(time),bc_Z(bc_Z)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV,class T2,class Z_INTERP,class T_AVERAGING,class U_INTERP,class T_FACE_LOOKUP>
ADVECTION_UPWIND<TV,T2,Z_INTERP,T_AVERAGING,U_INTERP,T_FACE_LOOKUP>::
~ADVECTION_UPWIND()
{
}
//#####################################################################
// Function Clamped_To_Array
//#####################################################################
template<class TV,class T2,class Z_INTERP,class T_AVERAGING,class U_INTERP,class T_FACE_LOOKUP>
T2 ADVECTION_UPWIND<TV,T2,Z_INTERP,T_AVERAGING,U_INTERP,T_FACE_LOOKUP>::
Clamped_To_Array(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,TV_INT>& Z,const TV& X) const
{
    T offset=grid.dX.Max();
    T phi_X=levelset.Phi(X);
    if(phi_X<-max_in)
    {
        Add_Debug_Particle(X,VECTOR<T,3>(0,1,0));
        return Z_INTERP().Clamped_To_Array(grid,Z,X);
    }
    TV N=levelset.Normal(X);
    T a=0,b=0;
    if(phi_X<=0){
        if(levelset.Phi(X+max_in*N)<=0)
        {
            Add_Debug_Particle(X,VECTOR<T,3>(0,1,1));
            return Z_INTERP().Clamped_To_Array(grid,Z,X);
        }
        b=max_in;}
    else{
        if(levelset.Phi(X-max_out*N)>0)
        {
            Add_Debug_Particle(X,VECTOR<T,3>(1,0,0));
            return T2();
        }
        a=-max_out;}

    struct LINE_PHI:public NONLINEAR_FUNCTION<T(T)>
    {
        const LEVELSET<TV>* levelset;
        TV X,dX;
        void Compute(const T t,T* ddf,T* df,T* f) const
        {
            TV Y=X+dX*t;
            if(f) *f=levelset->Phi(Y);
            if(df) *df=levelset->Normal(Y).Dot(dX);
            if(ddf) *ddf=(levelset->Hessian(Y)*dX).Dot(dX);
        }
    } line;
    line.levelset=&levelset;
    line.X=X;
    line.dX=N;

    ITERATIVE_SOLVER<T> iterative_solver;
    iterative_solver.tolerance=(T)1e-5;
    T c=iterative_solver.Bisection_Secant_Root(line,a,b);
    TV Y=X+N*c;
    TV V=U_INTERP().Clamped_To_Array_Face(grid,face_velocities,Y);
    TV NY=levelset.Normal(Y);
    if(V.Dot(NY)>0)
    {
        Add_Debug_Particle(X,VECTOR<T,3>(1,1,0));
        return Z_INTERP().Clamped_To_Array(grid,Z,X);
    }
    Add_Debug_Particle(Y,VECTOR<T,3>(.5,.5,.5));
    Add_Debug_Particle(Y-N*offset,VECTOR<T,3>(1,1,1));
    Add_Debug_Particle(Y-N*offset*2,VECTOR<T,3>(1,1,1));

    T2 A=bc_Z(Y+V*dt,time);
    T2 B=Z_INTERP().Clamped_To_Array(grid,Z,Y-N*offset);
    T2 C=Z_INTERP().Clamped_To_Array(grid,Z,Y-N*offset*2);
    T t=c/(2*offset);
    T2 D=A-t*(A*3+C-B*4)+2*t*t*(C+A-B*2);
//    D=A+c/offset*(B-A);
    Add_Debug_Particle(X,VECTOR<T,3>(1,0,1));
    return D;
}
//#####################################################################
// Function Update_Advection_Equation_Cell_Lookup
//#####################################################################
template<class TV,class T2,class Z_INTERP,class T_AVERAGING,class U_INTERP,class T_FACE_LOOKUP>
void ADVECTION_UPWIND<TV,T2,Z_INTERP,T_AVERAGING,U_INTERP,T_FACE_LOOKUP>::
Update_Advection_Equation_Cell_Lookup(const GRID<TV>& grid,ARRAY<T2,TV_INT>& Z,const ARRAY<T2,TV_INT>& Z_ghost,
    const T_FACE_LOOKUP& face_velocities,BOUNDARY<TV,T2>& boundary,const T dt,const T time,
    const ARRAY<T2,TV_INT>* Z_min_ghost,const ARRAY<T2,TV_INT>* Z_max_ghost,ARRAY<T2,TV_INT>* Z_min,ARRAY<T2,TV_INT>* Z_max)
{
    PHYSBAM_ASSERT(!Z_min && !Z_max);
    T_AVERAGING averaging;
    for(CELL_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()){
        TV_INT cell=iterator.Cell_Index();
        TV Y=iterator.Location()-dt*averaging.Face_To_Cell_Vector(grid,cell,face_velocities);
        Z(cell)=Clamped_To_Array(grid,Z_ghost,Y);}
}
template class ADVECTION_UPWIND<VECTOR<double,3>,SYMMETRIC_MATRIX<double,3>,
    QUADRATIC_INTERPOLATION_UNIFORM<VECTOR<double,3>,SYMMETRIC_MATRIX<double,3> > >;
template class ADVECTION_UPWIND<VECTOR<float,3>,SYMMETRIC_MATRIX<float,3>,
    QUADRATIC_INTERPOLATION_UNIFORM<VECTOR<float,3>,SYMMETRIC_MATRIX<float,3> > >;
template class ADVECTION_UPWIND<VECTOR<double,2>,SYMMETRIC_MATRIX<double,2>,
    QUADRATIC_INTERPOLATION_UNIFORM<VECTOR<double,2>,SYMMETRIC_MATRIX<double,2> > >;
template class ADVECTION_UPWIND<VECTOR<float,2>,SYMMETRIC_MATRIX<float,2>,
    QUADRATIC_INTERPOLATION_UNIFORM<VECTOR<float,2>,SYMMETRIC_MATRIX<float,2> > >;
}

