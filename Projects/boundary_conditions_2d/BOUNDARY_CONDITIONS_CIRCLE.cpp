#include <PhysBAM_Tools/Grids_Uniform/FACE_ITERATOR.h>
#include <climits>
#include <cstdio>
#include "BOUNDARY_CONDITIONS_CIRCLE.h"
#include "HEADER.h"

template<class TV> BOUNDARY_CONDITIONS_CIRCLE<TV>::
BOUNDARY_CONDITIONS_CIRCLE(GRID<TV>& grid_input,int bc_type)
    :BOUNDARY_CONDITIONS<TV>(grid_input),type(bc_type)
{
    bounding_box=base_domain=sphere.Bounding_Box();
    LOG::cout<<"HERE"<<std::endl;
    this->use_analytic_solution=true;
}

template<class TV> BOUNDARY_CONDITIONS_CIRCLE<TV>::
~BOUNDARY_CONDITIONS_CIRCLE()
{
    return;
    int n=0;
    T L1=0,Li=0;
    for(FACE_ITERATOR<TV> it(grid,2);it.Valid();it.Next()){
        TV X(grid.Face(it.Full_Index()));
        if(Theta(X)>grid.dX.Max()*2) continue;
        if(!Inside(X)) continue;
        TV Y;
        n++;
        T q=abs((*saved_u)(it.Full_Index())-Y(it.Axis()));
        TV dX=X-bounding_box.Center();
        printf("%.15g %.15g -> %.15g %.15g %.15g\n",(*saved_u)(it.Full_Index()),Y(it.Axis()),q,Theta(X),atan2(dX.y,dX.x));
        L1+=q;
        Li=max(Li,q);}

    LOG::cout<<"W@Z-L1 "<<grid.counts.x<<" "<<L1/n<<std::endl;
    LOG::cout<<"W@Z-Linf "<<grid.counts.x<<" "<<Li<<std::endl;
}

template<class TV> typename TV::SCALAR BOUNDARY_CONDITIONS_CIRCLE<TV>::
Theta(const TV& X) const
{
    return sphere.Signed_Distance(X);
}

template<class TV> bool BOUNDARY_CONDITIONS_CIRCLE<TV>::
Inside(const TV& X) const
{
    return sphere.Lazy_Inside(X);
}

template<class TV> int BOUNDARY_CONDITIONS_CIRCLE<TV>::
Boundary_Condition(const TV& X,const TV& Y,T& theta,TV& value,T time) const
{
    T py=sphere.Signed_Distance(Y),px=sphere.Signed_Distance(X);
    theta=px/(px-py);
    value=TV();
    if(type==source) PHYSBAM_NOT_IMPLEMENTED();
    return type;
}

template<class TV> void BOUNDARY_CONDITIONS_CIRCLE<TV>::
Initialize_Velocity_Field(ARRAY<T,FACE_INDEX<TV::m> >& u,T time) const
{
    saved_u=&u;
    for(FACE_ITERATOR<TV> it(grid,2);it.Valid();it.Next())
        u(it.Full_Index())=Analytic_Velocity(it.Location(),time)(it.Axis());
    Add_Initial_Error(u,time);
}

template<class TV> void BOUNDARY_CONDITIONS_CIRCLE<TV>::
Add_Initial_Error(ARRAY<T,FACE_INDEX<TV::m> >& u,T time) const
{
    for(FACE_ITERATOR<TV> it(grid,2);it.Valid();it.Next()){TV X=it.Location();
        TV V((X.x*X.x-X.x)*(X.y*X.y*X.y/3-X.y*X.y/2),(X.y*X.y-X.y)*(X.x*X.x*X.x/3-X.x*X.x/2));
        u(it.Full_Index())+=V(it.Axis());}
}

template<class TV> TV BOUNDARY_CONDITIONS_CIRCLE<TV>::
Analytic_Velocity(const TV& X,T time) const
{
    return X.Orthogonal_Vector()*0;
}

template struct BOUNDARY_CONDITIONS_CIRCLE<VECTOR<double,2> >;
