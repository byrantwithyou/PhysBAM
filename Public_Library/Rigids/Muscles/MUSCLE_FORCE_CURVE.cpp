//#####################################################################
// Copyright 2005, Ron Fedkiw, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Rigids/Muscles/MUSCLE_FORCE_CURVE.h>
using namespace PhysBAM;
template<class T> MUSCLE_FORCE_CURVE<T>::
MUSCLE_FORCE_CURVE()
    :passive_interpolation(&default_interpolation),active_interpolation(&default_interpolation),tendon_force_interpolation(&default_interpolation),velocity_interpolation(&default_interpolation)
{
}
template<class T> MUSCLE_FORCE_CURVE<T>::
~MUSCLE_FORCE_CURVE()
{
}
namespace{
//#####################################################################
// Function Compute_Inverse_Map_Helper
//#####################################################################
template<class T>
void Compute_Inverse_Map_Helper(const GRID<VECTOR<T,1> >& domain_grid,const ARRAY<T,VECTOR<int,1> >& function,const GRID<VECTOR<T,1> >& range_grid,ARRAY<T,VECTOR<int,1> >& inverse_function)
{
    int domain_i=1;
    T xmin=domain_grid.X(VECTOR<int,1>(0)).x,xmax=domain_grid.X(VECTOR<int,1>(domain_grid.counts.x)).x;
    for(int i=0;i<range_grid.counts.x;i++){
        T function_value=range_grid.X(VECTOR<int,1>(i)).x;
        while(domain_i<domain_grid.counts.x-1 && function(domain_i+1)<function_value) domain_i++;
        inverse_function(i)=clamp(domain_grid.X(VECTOR<int,1>(domain_i)).x+domain_grid.dX.x*(function_value-function(domain_i))/(function(domain_i+1)-function(domain_i)),xmin,xmax);}
}
//#####################################################################
// Function Compute_Inverse_Map_Helper
//#####################################################################
template<class T,class T2,int d>
void Compute_Inverse_Map_Helper(const GRID<VECTOR<T,d> >& domain_grid,const ARRAY<T2,VECTOR<int,1> >& function,const GRID<VECTOR<T,d> >& range_grid,ARRAY<T2,VECTOR<int,1> >& inverse_function)
{
    PHYSBAM_FATAL_ERROR();
}
}
template<class T> void MUSCLE_FORCE_CURVE<T>::
Initialize(const std::string& data_directory)
{
    Load_Data(data_directory+"/Muscle_Curves/passive_cubic",passive_force_grid,passive_force);
    Load_Data(data_directory+"/Muscle_Curves/active_cubic",active_force_grid,active_force);
    Load_Data(data_directory+"/Muscle_Curves/tendon_cubic",tendon_force_grid,tendon_force);
    velocity_grid.Initialize(VECTOR<int,1>(1001),RANGE<VECTOR<T,1> >(VECTOR<T,1>(-2),VECTOR<T,1>(2)));
    velocity_curve.Resize(velocity_grid.Domain_Indices());
    for(int v=0;v<1001;v++){
        T velocity=velocity_grid.X(VECTOR<int,1>(v)).x;
        velocity_curve(v)=(velocity<-1)?0:(T).54*atan((T)5.69*velocity+(T).51)+(T).745;}
    tendon_length_grid.Initialize(VECTOR<int,1>(1001),RANGE<VECTOR<T,1> >(VECTOR<T,1>(),VECTOR<T,1>(3.5)));tendon_length.Resize(tendon_length_grid.Domain_Indices());
    Compute_Inverse_Map_Helper(tendon_force_grid,tendon_force,tendon_length_grid,tendon_length);
    Compute_Slopes(passive_force_grid,passive_force,passive_force_slope_grid,passive_force_slope);
    Compute_Slopes(active_force_grid,active_force,active_force_slope_grid,active_force_slope);
    Compute_Slopes(tendon_force_grid,tendon_force,tendon_force_slope_grid,tendon_force_slope);
}
/*function above from "Model the leg cycling movement with neural oscillator", Dingguo Zhang; Kuanyi Zhu; Hang Zheng;Systems, Man and Cybernetics, 2004 IEEE International Conference on Volume 1,  10-13 Oct. 2004 Page(s):740 - 744 vol.1 */
template<class T> void MUSCLE_FORCE_CURVE<T>::
Load_Data(const std::string& prefix,GRID<VECTOR<T,1> >& grid,ARRAY<T,VECTOR<int,1> >& values)
{
    {std::istream* input=FILE_UTILITIES::Safe_Open_Input(prefix+".grid",false);
    *input>>grid.counts.x>>grid.domain.min_corner.x>>grid.domain.max_corner.x;
    grid.Initialize(grid.counts,grid.domain);
    values.Resize(grid.Domain_Indices());
    delete input;}
    std::istream* input(FILE_UTILITIES::Safe_Open_Input(prefix+".values",false));
    for(int i=0;i<values.array.m;i++) *input>>values.array(i);
    delete input;
}
template<class T> void MUSCLE_FORCE_CURVE<T>::
Compute_Slopes(const GRID<VECTOR<T,1> >& grid,const ARRAY<T,VECTOR<int,1> >& values,GRID<VECTOR<T,1> >& slope_grid,ARRAY<T,VECTOR<int,1> >& slopes)
{
    assert(!grid.Is_MAC_Grid());
    slope_grid=grid.Get_MAC_Grid();
    slopes.Resize(0,slope_grid.counts.x);
    for(int i=0;i<slopes.array.m;i++) slopes.array(i)=(values.array(i+1)-values.array(i))*grid.one_over_dX.x;
}
namespace PhysBAM{
template class MUSCLE_FORCE_CURVE<float>;
template class MUSCLE_FORCE_CURVE<double>;
}
