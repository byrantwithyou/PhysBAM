//#####################################################################
// Copyright 2002-2005, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EXTRAPOLATION
//#####################################################################
#ifndef __EXTRAPOLATION__
#define __EXTRAPOLATION__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Math_Tools/constants.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
namespace PhysBAM{

template<class T_GRID> struct GRID_ARRAYS_POLICY;
template<class T_GRID,class T2> class BOUNDARY_UNIFORM;

template<class T_GRID,class T2>
class EXTRAPOLATION:public NONCOPYABLE
{
    typedef typename T_GRID::SCALAR T;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_BASE T_ARRAYS_BASE;typedef typename T_ARRAYS_BASE::template REBIND<bool>::TYPE T_ARRAYS_BOOL_BASE;
    typedef typename T_GRID::INDEX T_INDEX;
public:
    T band_width; // band for extrapolation near the interface
    T isobaric_fix_width;  // band for the isobaric fix
    BOUNDARY_UNIFORM<T_GRID,T2>* boundary;
private:
    BOUNDARY_UNIFORM<T_GRID,T2>& boundary_default;

protected:
    EXTRAPOLATION()
        :boundary_default(*new BOUNDARY_UNIFORM<T_GRID,T2>)
    {boundary=&boundary_default;}

    ~EXTRAPOLATION()
    {delete &boundary_default;}

public:
    void Set_Custom_Boundary(BOUNDARY_UNIFORM<T_GRID,T2>& boundary_input)
    {boundary=&boundary_input;}

protected:
    static void Add_To_Heap(const T_ARRAYS_BASE& phi,ARRAY<T_INDEX>& heap,int& heap_length,T_ARRAYS_BOOL_BASE& close,const T_INDEX& index)
    {close(index)=true;heap(heap_length++)=index;Up_Heap(phi,heap,heap_length-1);}

    static T_INDEX Remove_Root_From_Heap(const T_ARRAYS_BASE& phi,ARRAY<T_INDEX>& heap,int& heap_length,T_ARRAYS_BOOL_BASE& close)
    {T_INDEX index=heap(0);close(index)=false;Down_Heap(phi,heap,heap_length);heap_length--;return index;}

    static void Up_Heap(const T_ARRAYS_BASE& phi,ARRAY<T_INDEX>& heap,int child)
    {int parent=(child-1)/2; // k=(2k+1)/2 and k=(2k+2)/2 from integer division
    while(parent >= 0 && phi(heap(child)) < phi(heap(parent))){
        exchange(heap(child),heap(parent)); // exchange child & parent
        child=parent;parent=(child-1)/2;}} // move up one

    static void Down_Heap(const T_ARRAYS_BASE& phi,ARRAY<T_INDEX>& heap,const int heap_length)
    {int parent=0,child_left=1,child_right=2; // initialize
    while(heap_length > child_right){
        if(phi(heap(child_left)) <= phi(heap(child_right))){
            heap(parent)=heap(child_left); // update index
            parent=child_left;}  // move down one
        else{
            heap(parent)=heap(child_right); // update index
            parent=child_right;} // move down one
        child_left=2*parent+1;child_right=2*parent+2;} // find new children
    // fill the hole with the last element
    if(parent != heap_length-1){
        heap(parent)=heap(heap_length-1); // update index
        Up_Heap(phi,heap,parent);}}

    static T2 Solve_Close_Point(const bool no_x,const T2& value_x,const T phix_dx,const bool no_y,const T2& value_y,const T phiy_dy,const T dy_squared_over_dx_squared,const T tolerance)
    {assert(!no_x || !no_y);
    if(no_x) return value_y;
    if(no_y) return value_x;
    // candidates exist in both directions
    T a=phix_dx*dy_squared_over_dx_squared,b=phiy_dy,denominator=a+b,fraction=.5;
    if(denominator > tolerance) fraction=clamp(a/denominator,(T)0,(T)1);
    return fraction*value_x+(1-fraction)*value_y;}

    static T2 Solve_Close_Point(const bool no_x,const T2& value_x,const T phix_dx,const bool no_y,const T2& value_y,const T phiy_dy,const bool no_z,const T2& value_z,const T phiz_dz,
        const T dz_squared_over_dy_squared,const T dz_squared_over_dx_squared,const T dy_squared_over_dx_squared,const T tolerance)
    {assert(!no_x || !no_y || !no_z);
    if(no_x) return Solve_Close_Point(no_y,value_y,phiy_dy,no_z,value_z,phiz_dz,dz_squared_over_dy_squared,tolerance);
    if(no_y) return Solve_Close_Point(no_x,value_x,phix_dx,no_z,value_z,phiz_dz,dz_squared_over_dx_squared,tolerance);
    if(no_z) return Solve_Close_Point(no_x,value_x,phix_dx,no_y,value_y,phiy_dy,dy_squared_over_dx_squared,tolerance);
    // candidates exist in all three directions
    T a=phix_dx*dz_squared_over_dx_squared,b=phiy_dy*dz_squared_over_dy_squared,c=phiz_dz,denominator=a+b+c,fraction_1=(T)one_third,fraction_2=(T)one_third;
    if(denominator > tolerance){fraction_1=clamp(a/denominator,(T)0,(T)1);fraction_2=clamp(b/denominator,(T)0,1-fraction_1);}
    return fraction_1*value_x+fraction_2*value_y+(1-fraction_1-fraction_2)*value_z;}

//#####################################################################
};
}
#endif
