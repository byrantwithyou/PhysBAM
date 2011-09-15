#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Tools/Matrices/MATRIX.h>
#include "BOUNDARY_CONDITIONS_SEGMENT.h"
#include "PARAMETERS_COMMON.h"

template<class TV> BOUNDARY_CONDITIONS_SEGMENT<TV>::
BOUNDARY_CONDITIONS_SEGMENT(GRID<TV>& grid_input,T a,T b)
    :BOUNDARY_CONDITIONS<TV>(grid_input)
{
    PHYSBAM_ASSERT(a<=0 && b>=1);
    base_domain=RANGE<TV>::Unit_Box();
    bounding_box=RANGE<TV>(TV(a),TV(b));
    use_analytic_solution=true;
}

template<class TV> BOUNDARY_CONDITIONS_SEGMENT<TV>::
~BOUNDARY_CONDITIONS_SEGMENT()
{
}

template<class TV> typename TV::SCALAR BOUNDARY_CONDITIONS_SEGMENT<TV>::
Theta(const TV& X) const
{
    return (bounding_box.Center()-X).Magnitude()-bounding_box.Edge_Lengths().x/2;
}

template<class TV> int BOUNDARY_CONDITIONS_SEGMENT<TV>::
Boundary_Condition(const TV& X,const TV& Y,T& theta,TV& value,T time) const
{
    T py=Theta(Y),px=Theta(X);
    theta=px/(px-py);
    TV Z=X+(Y-X)*theta;
    value=Analytic_Velocity(Z,time);
    return BOUNDARY_CONDITIONS<TV>::noslip;
}

template<class TV> void BOUNDARY_CONDITIONS_SEGMENT<TV>::
Initialize_Velocity_Field(ARRAY<T,FACE_INDEX<TV::m> >& u,T time) const
{
    for(UNIFORM_GRID_ITERATOR_FACE<TV> it(grid,2);it.Valid();it.Next())
        u(it.Full_Index())=Analytic_Velocity(it.Location(),time)(it.Axis());
}

template<class TV> TV BOUNDARY_CONDITIONS_SEGMENT<TV>::
Analytic_Velocity(const TV& X,T time) const
{
    return sin((T)2*X)*exp(-4*nu*time);
}

template<class TV> void BOUNDARY_CONDITIONS_SEGMENT<TV>::
Update_Parameters(PARAMETERS_COMMON<T>& param)
{
    nu=param.mu/param.rho;
}

template class BOUNDARY_CONDITIONS_SEGMENT<VECTOR<double,1> >;
