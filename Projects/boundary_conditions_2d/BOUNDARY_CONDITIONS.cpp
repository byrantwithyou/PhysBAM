#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_1D.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_2D.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_3D.h>
#include "BOUNDARY_CONDITIONS.h"

template<class TV> BOUNDARY_CONDITIONS<TV>::
BOUNDARY_CONDITIONS(GRID<TV>& grid_input)
:grid(grid_input),check_leaks(false),use_analytic_solution(false),phi(0)
{
}
template<class TV> BOUNDARY_CONDITIONS<TV>::
~BOUNDARY_CONDITIONS()
{
}
template<class TV> TV BOUNDARY_CONDITIONS<TV>::
Gradient(const TV& X) const
{
    T dx=(T)1e-4;
    TV grad;
    for(int a=0;a<TV::m;a++){
        TV Y=X;Y(a)-=dx;
        T l=Theta(Y);
        Y(a)+=2*dx;
        T r=Theta(Y);
        grad(a)=(r-l)/(2*dx);}
    return grad;
}
template<class TV> bool BOUNDARY_CONDITIONS<TV>::
Inside(const TV& X) const
{
    return Theta(X)<=0;
}
template<class TV> void BOUNDARY_CONDITIONS<TV>::
Initialize_Phi(int ghost)
{
    phi_array.Resize(grid.Domain_Indices(ghost));
    for(UNIFORM_GRID_ITERATOR_CELL<TV> it(grid,ghost);it.Valid();it.Next())
        phi_array(it.index)=Theta(it.Location());
    if(!phi) phi=new typename LEVELSET_POLICY<GRID<TV> >::LEVELSET(grid,phi_array);
}

template<class TV> void BOUNDARY_CONDITIONS<TV>::
Compute_Error(ARRAY<T,FACE_INDEX<TV::m> >& u,T time) const
{
    if(!use_analytic_solution) return;
    int n=0;
    T L1=0,Li=0;
    for(UNIFORM_GRID_ITERATOR_FACE<TV> it(grid);it.Valid();it.Next()){
        TV X(grid.Axis_X_Face(it.Full_Index()));
        if(!Inside(X)) continue;
        TV p=Analytic_Velocity(X,time);
        n++;
        T q=abs(u(it.Full_Index())-p(it.Axis()));
        L1+=q;
        Li=max(Li,q);}

    LOG::cout<<"W@Z-L1 "<<grid.counts.x<<" "<<L1/n<<std::endl;
    LOG::cout<<"W@Z-Linf "<<grid.counts.x<<" "<<Li<<std::endl;
}

template<class TV> void BOUNDARY_CONDITIONS<TV>::
Update_Parameters(PARAMETERS_COMMON<T>& param)
{
}

template struct BOUNDARY_CONDITIONS<VECTOR<double,1> >;
template struct BOUNDARY_CONDITIONS<VECTOR<double,2> >;
