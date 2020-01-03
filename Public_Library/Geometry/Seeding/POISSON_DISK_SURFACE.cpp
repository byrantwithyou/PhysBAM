//#####################################################################
// Copyright 2020, Craig Schroeder, Yunxin Sun.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class POISSON_DISK_SURFACE
//#####################################################################
#include <Grid_Tools/Grids/GRID.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>
#include <Geometry/Seeding/POISSON_DISK_SURFACE.h>
using namespace PhysBAM;
//#####################################################################
// Function Sample
//#####################################################################
template<class TV> void POISSON_DISK_SURFACE<TV>::
Sample(RANDOM_NUMBERS<T>& random,IMPLICIT_OBJECT<TV>* io,ARRAY<TV>& X)
{
    PHYSBAM_ASSERT(h && min_distance);
    io->Update_Box();
    GRID<TV> grid=GRID<TV>::Create_Grid_Given_Cell_Size(io->box,h,true,0);

    ARRAY<int,TV_INT> grid_array(grid.Cell_Indices(ghost),use_init,-1);
    if(X.m){
        for(int i=0;i<X.m;i++)
            grid_array(grid.Cell(X(i)))=i;}
    else{
        TV first_point;
        random.Fill_Uniform(first_point,io->box);
        first_point=io->Closest_Point_On_Boundary(first_point);
        grid_array(grid.Cell(first_point))=X.Append(first_point);}
    ARRAY<int> active(IDENTITY_ARRAY<>(X.m));
    while(active.m){
        int random_index=random.Get_Uniform_Integer(0,active.m-1);
        TV point=X(active(random_index));
        bool found_at_least_one=false;
        for(int i=0;i<max_attemps;i++){
            TV new_point=Generate_Random_Point_Around_Annulus(random,point);
            if(!Adjust_To_Surface(io,point,new_point)) continue;
            TV_INT cell_index=grid.Cell(new_point);
            int& ci=grid_array(cell_index);
            if(ci!=-1) continue;
            if(!Check_Distance(grid,grid_array,new_point,X)) continue;
            found_at_least_one=true;
            int index=X.Append(new_point);
            active.Append(index);
            ci=index;
            break;}
        if(!found_at_least_one) active.Remove_Index_Lazy(random_index);}
}
//#####################################################################
// Function Generate_Random_Point_Around_Annulus
//#####################################################################
template<class TV> TV POISSON_DISK_SURFACE<TV>::
Generate_Random_Point_Around_Annulus(RANDOM_NUMBERS<T>& random,const TV& center) const
{
    while(1){
        TV v;
        random.Fill_Uniform(v,-1,1);
        T mag2=v.Magnitude_Squared();
        if(mag2>=0.25 && mag2<=1) return v*min_distance*2+center;}
    return TV();
}
//#####################################################################
// Function Generate_Random_Point_Around_Annulus
//#####################################################################
template<class TV> bool POISSON_DISK_SURFACE<TV>::
Adjust_To_Surface(IMPLICIT_OBJECT<TV>* io,const TV& seed,TV& pt) const
{
    for(int i=0;i<max_proj_attemps;i++){
        pt=io->Closest_Point_On_Boundary(pt);
        TV DX=pt-seed;
        T d=DX.Normalize();
        if(d>=min_distance && d<=2*min_distance) return true;
        pt=DX*(T)1.5*min_distance+seed;}
    return false;
}
//#####################################################################
// Function Check_Distance
//#####################################################################
template<class TV> bool POISSON_DISK_SURFACE<TV>::
Check_Distance(const GRID<TV>& grid,ARRAY<int,TV_INT>& grid_array,const TV& point,ARRAY<TV>& X) const
{
    TV_INT cell=grid.Cell(point);
    RANGE<TV_INT> candidate_range(cell-2,cell+3);
    for(RANGE_ITERATOR<TV::m> it(candidate_range);it.Valid();it.Next()){
        int i=grid_array(it.index);
        if(i>=0 && (point-X(i)).Magnitude_Squared()<sqr(min_distance)) return false;}
    return true;
}
//#####################################################################
// Function Set_Distance_By_Volume
//#####################################################################
template<class TV> void POISSON_DISK_SURFACE<TV>::
Set_Distance_By_Volume(T volume_per_sample)
{
    if(TV::m==1) min_distance=volume_per_sample*((T)2/3);
    else if(TV::m==2) min_distance=sqrt(volume_per_sample*((T)2/3));
    else min_distance=cbrt(volume_per_sample*((T)13/18));
    h=min_distance/sqrt((T)(TV::m));
}
//#####################################################################
// Function Set_Distance
//#####################################################################
template<class TV> void POISSON_DISK_SURFACE<TV>::
Set_Distance(T distance)
{
    min_distance=distance;
    h=distance/sqrt((T)(TV::m));
}
namespace PhysBAM{
template class POISSON_DISK_SURFACE<VECTOR<float,1> >;
template class POISSON_DISK_SURFACE<VECTOR<float,2> >;
template class POISSON_DISK_SURFACE<VECTOR<float,3> >;
template class POISSON_DISK_SURFACE<VECTOR<double,1> >;
template class POISSON_DISK_SURFACE<VECTOR<double,2> >;
template class POISSON_DISK_SURFACE<VECTOR<double,3> >;
}
