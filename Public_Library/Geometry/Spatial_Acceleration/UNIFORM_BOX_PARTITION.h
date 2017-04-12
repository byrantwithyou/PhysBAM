//#####################################################################
// Copyright 2004-2005, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class UNIFORM_BOX_PARTITION
//#####################################################################
#ifndef __UNIFORM_BOX_PARTITION__
#define __UNIFORM_BOX_PARTITION__

#include <Core/Arrays/ARRAY.h>
#include <Core/Data_Structures/PAIR.h>
#include <Core/Utilities/NONCOPYABLE.h>
#include <Grid_Tools/Grids/GRID.h>
#include <Geometry/Basic_Geometry/RAY.h>
#include <Geometry/Intersections/RAY_BOX_INTERSECTION.h>
namespace PhysBAM{

template<class T,class DATA_T>
class UNIFORM_BOX_PARTITION:public NONCOPYABLE
{
    typedef VECTOR<T,3> TV;typedef VECTOR<int,3> TV_INT;
public:
    GRID<TV> grid;
    ARRAY<ARRAY<DATA_T>*,VECTOR<int,3> > cells;
    RANGE<TV> bounding_box;
    bool initialized;
    
    UNIFORM_BOX_PARTITION()
        :initialized(false)
    {}

    ~UNIFORM_BOX_PARTITION()
    {for(int i=0;i<cells.array.m;i++) delete cells.array(i);}

    void Initialize(ARRAY<PAIR<RANGE<TV>,DATA_T> >& boxes_input,const T thickness_over_two=1e-6)
    {if(boxes_input.m==0){grid.Initialize(TV_INT()+2,RANGE<TV>::Unit_Box());cells.Resize(grid.Domain_Indices());return;}
    bounding_box=boxes_input(0).x;for(int k=1;k<boxes_input.m;k++) bounding_box.Enlarge_To_Include_Box(boxes_input(k).x.Thickened(thickness_over_two));
    bounding_box=bounding_box.Thickened(thickness_over_two);
    VECTOR<T,3> lengths=bounding_box.Edge_Lengths();
    VECTOR<int,3> dimensions;
    int max_axis=lengths.Dominant_Axis();T max_dimension=((T)3*pow((T)boxes_input.m,((T)1/3)));
    dimensions[max_axis]=(int)max_dimension; 
    int other_axis_1=(max_axis+1)%3;int other_axis_2=(other_axis_1+1)%3; // rotate to get other axis indices
    dimensions[other_axis_1]=(int)(max_dimension*T(lengths[other_axis_1])/T(lengths[max_axis]));
    dimensions[other_axis_2]=(int)(max_dimension*T(lengths[other_axis_2])/T(lengths[max_axis]));
    dimensions=clamp_min(dimensions,VECTOR<int,3>(1,1,1));
    grid.Initialize(dimensions+1,bounding_box);
    if(initialized) for(int i=0;i<cells.array.m;i++) delete cells.array(i);
    cells.Resize(grid.Get_MAC_Grid().Domain_Indices(),false,false);cells.Fill(0);initialized=true;
    for(int k=0;k<boxes_input.m;k++){
        RANGE<TV_INT> domain(grid.Cell(boxes_input(k).x.min_corner),grid.Cell(boxes_input(k).x.max_corner));
        domain.max_corner=TV_INT::Componentwise_Min(domain.max_corner,grid.numbers_of_cells);
        for(RANGE_ITERATOR<TV::m> it(domain);it.Valid();it.Next()){
            if(!cells(it.index))cells(it.index)=new ARRAY<DATA_T>;
            cells(it.index)->Append(boxes_input(k).y);}}}

//#####################################################################
// Function Map_Intersection
//#####################################################################
template<class HELPER_T>
bool Map_Intersection(RAY<VECTOR<T,3> >& ray,HELPER_T pointer) const
{
    T t_start=0;
    if(bounding_box.Lazy_Outside(ray.endpoint)){RAY<VECTOR<T,3> > ray_temp=ray;if(INTERSECTION::Intersects(ray_temp,bounding_box)) t_start=ray_temp.t_max;else return false;}
    VECTOR<T,3> point=ray.Point(t_start);
    VECTOR<int,3> index=grid.Clamped_Index_End_Minus_One(point);
    VECTOR<T,3> cross,dt;VECTOR<int,3> step,end;T parallel_tolerance=(T)1e-6*grid.dX.Min();
    TV X=grid.X(index),X0=grid.X(index+1);
    if(abs(ray.direction.x) < parallel_tolerance){cross.x=FLT_MAX;dt.x=0;end.x=0;}
    else{
        T one_over_direction_x=1/ray.direction.x;
        if(ray.direction.x > 0){cross.x=t_start+(X0.x-point.x)*one_over_direction_x;dt.x=grid.dX.x*one_over_direction_x;step.x=1;end.x=grid.counts.x;}
        else{cross.x=t_start+(X.x-point.x)*one_over_direction_x;dt.x=-grid.dX.x*one_over_direction_x;step.x=-1;end.x=0;}}
    if(abs(ray.direction.y) < parallel_tolerance){cross.y=FLT_MAX;dt.y=0;end.y=0;}
    else{
        T one_over_direction_y=1/ray.direction.y;
        if(ray.direction.y > 0){cross.y=t_start+(X0.y-point.y)*one_over_direction_y;dt.y=grid.dX.y*one_over_direction_y;step.y=1;end.y=grid.counts.y;}
        else{cross.y=t_start+(X.y-point.y)*one_over_direction_y;dt.y=-grid.dX.y*one_over_direction_y;step.y=-1;end.y=0;}}
    if(abs(ray.direction.z) < parallel_tolerance){cross.z=FLT_MAX;dt.z=0;end.z=0;}
    else{
        T one_over_direction_z=1/ray.direction.z;
        if(ray.direction.z > 0){cross.z=t_start+(X0.z-point.z)*one_over_direction_z;dt.z=grid.dX.z*one_over_direction_z;step.z=1;end.z=grid.counts.z;}
        else{cross.z=t_start+(X.z-point.z)*one_over_direction_z;dt.z=-grid.dX.z*one_over_direction_z;step.z=-1;end.z=0;}}
    bool intersection=false;
    for(;;){
        int compare_bits=((cross[0]<cross[1])<<2)+((cross[0]<cross[2])<<1)+((cross[1]<cross[2]));
        const int compare_bits_to_axis[8]={3,2,3,2,3,3,1,1};int axis=compare_bits_to_axis[compare_bits];
        if(cells(index)&&pointer->Callback(ray,*cells(index),cross[axis]))return true;
        if(!ray.semi_infinite&&ray.t_max<cross[axis])break;
        index[axis]+=step[axis];
        if(index[axis]==end[axis])break;
        cross[axis]+=dt[axis];}
    return intersection;
}
//#####################################################################
};
}
#endif
