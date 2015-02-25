//#####################################################################
// Copyright 2015, Andre Pradhana, Chenfanfu Jiang.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class POISSON_DISK
//#####################################################################
#include <Tools/Grids_Uniform/GRID.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>
#include <Geometry/Seeding/POISSON_DISK.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> POISSON_DISK<TV>::
POISSON_DISK(const typename TV::SCALAR min_distance,const int max_attemps)
    :min_distance(min_distance),max_attemps(max_attemps),ghost(6),h(min_distance/sqrt((T)(TV::m)))
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> POISSON_DISK<TV>::
~POISSON_DISK()
{
}
//#####################################################################
// Function Sample
//#####################################################################
template<class TV> void POISSON_DISK<TV>::
Sample(IMPLICIT_OBJECT<TV>* object,ARRAY<TV>& XXX)
{
    object->Update_Box();
    RANGE<TV> bounding_box=object->box.Thickened(min_distance*2);
    TV_INT cell_counts=TV_INT(bounding_box.Edge_Lengths()/h)+1;
    GRID<TV> grid(cell_counts+1,RANGE<TV>(bounding_box.min_corner,bounding_box.min_corner+(TV)cell_counts*h));
    ARRAY<int,TV_INT> grid_array(grid.Cell_Indices(ghost),true,-1);
    PHYSBAM_ASSERT(grid.domain.Contains(object->box));
    ARRAY<TV> XX;
    ARRAY<int> active;
    TV first_point=grid.domain.Center();
    grid_array(grid.Cell(first_point,ghost))=0;
    XX.Append(first_point);
    active.Append(0);
    while(active.m){
        int random_index=ran.Get_Uniform_Integer(0,active.m-1);
        TV point=XX(active(random_index));
        bool found_at_least_one=false;
        for(int i=0;i<max_attemps;i++){
            TV new_point=Generate_Random_Point_Around_Annulus(ran,point);
            if(Check_Distance(grid,grid_array,new_point,XX) && grid.domain.Lazy_Inside(new_point)){
                found_at_least_one=true;
                int index=XX.Append(new_point);
                active.Append(index);
                grid_array(grid.Cell(new_point,ghost))=index;}}
        if(!found_at_least_one) active.Remove_Index_Lazy(random_index);}
    for(int i=0;i<XX.m;i++) if(object->Extended_Phi(XX(i))<=0) XXX.Append(XX(i));
}
//#####################################################################
// Function Generate_Random_Point_Around_Annulus
//#####################################################################
template<class TV> TV POISSON_DISK<TV>::
Generate_Random_Point_Around_Annulus(RANDOM_NUMBERS<T>& random,TV& center) const
{
    while(1){
        TV v=random.Get_Uniform_Vector(RANGE<TV>::Centered_Box());
        T mag2=v.Magnitude_Squared();
        if(mag2>=0.25 && mag2<=1) return v*min_distance*2+center;}
    return TV();
}
//#####################################################################
// Function Check_Distance
//#####################################################################
template<class TV> bool POISSON_DISK<TV>::
Check_Distance(const GRID<TV>& grid,ARRAY<int,TV_INT>& grid_array,const TV& point,ARRAY<TV>& X) const
{
    TV_INT cell=grid.Cell(point,ghost);
    RANGE<TV_INT> candidate_range(cell-1,cell+2);
    for(RANGE_ITERATOR<TV::m> it(candidate_range);it.Valid();it.Next()){
        if(grid_array(it.index)==-1) continue;
        if((point-X(grid_array(it.index))).Magnitude_Squared()<sqr(min_distance)) return false;}
    return true;
}
namespace PhysBAM{
template class POISSON_DISK<VECTOR<float,2> >;
template class POISSON_DISK<VECTOR<float,3> >;
template class POISSON_DISK<VECTOR<double,2> >;
template class POISSON_DISK<VECTOR<double,3> >;
}
