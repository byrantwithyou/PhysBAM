#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Tools/Matrices/MATRIX.h>
#include "BOUNDARY_CONDITIONS_GIBOU_ANALYTIC_ONE.h"

template<class TV> BOUNDARY_CONDITIONS_GIBOU_ANALYTIC_ONE<TV>::
BOUNDARY_CONDITIONS_GIBOU_ANALYTIC_ONE(GRID<TV>& grid_input)
    :BOUNDARY_CONDITIONS<TV>(grid_input)
{
    bounding_box=base_domain=RANGE<TV>::Unit_Box()*(T)pi;
    use_analytic_solution=true;
}

template<class TV> BOUNDARY_CONDITIONS_GIBOU_ANALYTIC_ONE<TV>::
~BOUNDARY_CONDITIONS_GIBOU_ANALYTIC_ONE()
{
}

template<class TV> typename TV::SCALAR BOUNDARY_CONDITIONS_GIBOU_ANALYTIC_ONE<TV>::
Theta(const TV& X) const
{
    TV w=(X-pi/2).Normalized()+pi/2;
    VECTOR<T,3> z(w.x,w.y,0);
    T k=(T).2;
    for(int i=1;i<100;i++){
        T cx=cos(z.x),cy=cos(z.y),sx=sin(z.x),sy=sin(z.y);
        VECTOR<T,3> G(-2*X.x+2*z.x+z.z*cx*sy,-2*X.y+2*z.y+z.z*sx*cy,-k+sx*sy);
        MATRIX<T,3> H(2-z.z*sx*sy,z.z*cx*cy,cx*sy,z.z*cx*cy,2-z.z*sx*sy,sx*cy,cx*sy,sx*cy,0);
        z-=H.Solve_Linear_System(G);
        if(G.Magnitude_Squared()<1e-25) break;}
    T sign=((T).2-sin(X.x)*sin(X.y))>0?1:-1;
    if(!bounding_box.Lazy_Inside(X)) sign=1;
    return (z.Remove_Index(2)-X).Magnitude()*sign;
}

template<class TV> int BOUNDARY_CONDITIONS_GIBOU_ANALYTIC_ONE<TV>::
Boundary_Condition(const TV& X,const TV& Y,T& theta,TV& value,T time) const
{
    T py=Theta(Y),px=Theta(X),pw=Theta(X*2-Y);
    T d=(py-pw)/2,dd=pw+py-2*px,den=d+sqrt(d*d-2*dd*px);
    theta=-2*px/den;
    TV Z=X+(Y-X)*theta;
    value=Analytic_Velocity(Z,time);
    return BOUNDARY_CONDITIONS<TV>::noslip;
}

template<class TV> void BOUNDARY_CONDITIONS_GIBOU_ANALYTIC_ONE<TV>::
Initialize_Velocity_Field(ARRAY<T,FACE_INDEX<TV::m> >& u,T time) const
{
    for(UNIFORM_GRID_ITERATOR_FACE<TV> it(grid,2);it.Valid();it.Next())
        u(it.Full_Index())=Analytic_Velocity(it.Location(),time)(it.Axis());
    Add_Initial_Error(u,time);
}

template<class TV> void BOUNDARY_CONDITIONS_GIBOU_ANALYTIC_ONE<TV>::
Add_Initial_Error(ARRAY<T,FACE_INDEX<TV::m> >& u,T time) const
{
    for(UNIFORM_GRID_ITERATOR_FACE<TV> it(grid,2);it.Valid();it.Next()){TV X=it.Location();
        TV V((X.x*X.x-X.x)*(X.y*X.y*X.y/3-X.y*X.y/2),(X.y*X.y-X.y)*(X.x*X.x*X.x/3-X.x*X.x/2));
        u(it.Full_Index())+=V(it.Axis());}
}

template<class TV> TV BOUNDARY_CONDITIONS_GIBOU_ANALYTIC_ONE<TV>::
Analytic_Velocity(const TV& X,T time) const
{
    return TV(sin(X.x)*cos(X.y),-cos(X.x)*sin(X.y))*cos(time);
}

template<class TV> void BOUNDARY_CONDITIONS_GIBOU_ANALYTIC_ONE<TV>::
Update_Parameters(PARAMETERS_COMMON<T>& param)
{
    param.mu=tan(param.time)*param.rho/2;
}

template struct BOUNDARY_CONDITIONS_GIBOU_ANALYTIC_ONE<VECTOR<double,2> >;
