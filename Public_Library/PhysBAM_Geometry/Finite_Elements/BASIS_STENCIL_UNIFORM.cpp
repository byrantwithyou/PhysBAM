//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Symbolics/STATIC_POLYNOMIAL.h>
#include <PhysBAM_Geometry/Finite_Elements/BASIS_STENCIL_UNIFORM.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV,int d> BASIS_STENCIL_UNIFORM<TV,d>::
BASIS_STENCIL_UNIFORM(TV dx)
    :dX(dx)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV,int d> BASIS_STENCIL_UNIFORM<TV,d>::
~BASIS_STENCIL_UNIFORM()
{
}
//#####################################################################
// Function Add_Symmetric_Entry
//#####################################################################
template<class TV,int d> void BASIS_STENCIL_UNIFORM<TV,d>::
Add_Symmetric_Entry(const ENTRY& e, int mask) // 1=x, 2=y, 4=z
{
    int m=stencils.m;
    stencils.Append(e);
    for(int i=0;i<TV::m;i++)
        if(mask&(1<<i))
            for(int j=m,n=stencils.m;j<n;j++){
                ENTRY f=stencils(j);
                f.region.min_corner(i)=-stencils(j).region.max_corner(i);
                f.region.max_corner(i)=-stencils(j).region.min_corner(i);
                f.polynomial.Negate_Variable(i);
                stencils.Append(f);}
}
//#####################################################################
// Function Exchange_Axes
//#####################################################################
template<class TV,int d> void BASIS_STENCIL_UNIFORM<TV,d>::
Exchange_Axes(int m,int n)
{
    if(m==n) return;
    for(int i=0;i<stencils.m;i++){
        exchange(stencils(i).region.min_corner(m),stencils(i).region.min_corner(n));
        exchange(stencils(i).region.max_corner(m),stencils(i).region.max_corner(n));
        stencils(i).polynomial.Exchange_Variables(m,n);}
}
//#####################################################################
// Function Scale_Axes
//#####################################################################
template<class TV,int d> void BASIS_STENCIL_UNIFORM<TV,d>::
Scale_Axes(TV scale)
{
    for(int i=0;i<stencils.m;i++)
        stencils(i).polynomial.Scale(scale);
}
//#####################################################################
// Function Dice_Stencil
//#####################################################################
template<class TV,int d> void BASIS_STENCIL_UNIFORM<TV,d>::
Dice_Stencil()
{
    const VECTOR<TV_INT,(1<<TV::m)>& counts=GRID<TV>::Binary_Counts(TV_INT());

    diced.Remove_All();
    const int big_shift=1<<30;
    const int half_big_shift=1<<29;
    RANGE<TV_INT> cell_box=RANGE<TV_INT>::Centered_Box();
    for(int i=0;i<stencils.m;i++){
        RANGE<TV_INT> offset_range=stencils(i).region.Translated(center_offset);
        RANGE<TV_INT> ra((big_shift+1-offset_range.max_corner)/2-half_big_shift,(big_shift+2-offset_range.min_corner)/2-half_big_shift);
        for(RANGE_ITERATOR<TV::m> it(ra);it.Valid();it.Next()){
            RANGE<TV_INT> cut_range=RANGE<TV_INT>::Intersect(offset_range.Translated(2*it.index),cell_box);
            DICED e={it.index};
            for(int b=0;b<(1<<TV::m);b++) if(cut_range.Lazy_Inside_Half_Open(counts(b)-1)) e.subcell|=1<<b;
            e.polynomial=stencils(i).polynomial;
            e.polynomial.Shift(-dX*(TV(it.index)+TV(center_offset)/2));
            diced.Append(e);}}
}
//#####################################################################
// Function Set_Constant_Stencil
//#####################################################################
template<class TV,int d> void BASIS_STENCIL_UNIFORM<TV,d>::
Set_Constant_Stencil()
{
    ENTRY e;
    e.region=RANGE<TV_INT>::Centered_Box();
    e.polynomial.Set_Term(TV_INT(),1);
    e.polynomial.Scale((T)1/dX);
    stencils.Append(e);
}
//#####################################################################
// Function Set_Multilinear_Stencil
//#####################################################################
template<class TV,int d> void BASIS_STENCIL_UNIFORM<TV,d>::
Set_Multilinear_Stencil()
{
    ENTRY e;
    e.region.min_corner=TV_INT();
    e.region.max_corner=TV_INT()+2;
    e.polynomial.Set_Term(TV_INT()+1,TV::m%2==0?1:-1);
    e.polynomial.Shift(TV()-1);
    e.polynomial.Scale((T)1/dX);
    Add_Symmetric_Entry(e);
}
//#####################################################################
// Function Differentiate
//#####################################################################
template<class TV,int d> void BASIS_STENCIL_UNIFORM<TV,d>::
Differentiate(int v)
{
    int k=0;
    for(int i=0;i<stencils.m;i++){
        stencils(i).polynomial=stencils(i).polynomial.Differentiate(v);
        stencils(i).polynomial.Compress_Size();
        if(stencils(i).polynomial.size!=TV_INT())
            stencils(k++)=stencils(i);}
    stencils.Resize(k);
}
//#####################################################################
// Function Print
//#####################################################################
template<class TV,int d> void BASIS_STENCIL_UNIFORM<TV,d>::
Print() const
{
    LOG::cout<<"Stencils: "<<center_offset<<std::endl;
    for(int i=0;i<stencils.m;i++)
        LOG::cout<<"  "<<stencils(i).region<<"  "<<stencils(i).polynomial<<std::endl;

    LOG::cout<<"Diced:"<<std::endl;
    for(int i=0;i<diced.m;i++)
        LOG::cout<<"  "<<diced(i).index_offset<<"  "<<diced(i).subcell<<"  "<<diced(i).polynomial<<std::endl;
}
//#####################################################################
// Function Padding
//#####################################################################
template<class TV,int d> int BASIS_STENCIL_UNIFORM<TV,d>::
Padding() const
{
    RANGE<TV_INT> support;
    for(int i=0;i<stencils.m;i++) support.Enlarge_To_Include_Box(stencils(i).region);
    return max(support.min_corner.Max_Abs(),support.max_corner.Max_Abs())/2;
}
//#####################################################################
// Function Overlap_Padding
//#####################################################################
template<class TV,int d> template<int d1> int BASIS_STENCIL_UNIFORM<TV,d>::
Overlap_Padding(const BASIS_STENCIL_UNIFORM<TV,d1>& s1) const
{
    int padding=0;
    for(int k=0;k<diced.m;k++)
        for(int l=0;l<s1.diced.m;l++)
            if(diced(k).subcell&s1.diced(l).subcell)
                padding=max((s1.diced(l).index_offset-diced(k).index_offset).Max_Abs(),padding);
    return padding;
}
template struct BASIS_STENCIL_UNIFORM<VECTOR<float,1>,0>;
template struct BASIS_STENCIL_UNIFORM<VECTOR<float,1>,1>;
template int BASIS_STENCIL_UNIFORM<VECTOR<float,1>,0>::Overlap_Padding<0>(BASIS_STENCIL_UNIFORM<VECTOR<float,1>,0> const&) const;
template int BASIS_STENCIL_UNIFORM<VECTOR<float,1>,0>::Overlap_Padding<1>(BASIS_STENCIL_UNIFORM<VECTOR<float,1>,1> const&) const;
template int BASIS_STENCIL_UNIFORM<VECTOR<float,1>,1>::Overlap_Padding<0>(BASIS_STENCIL_UNIFORM<VECTOR<float,1>,0> const&) const;
template int BASIS_STENCIL_UNIFORM<VECTOR<float,1>,1>::Overlap_Padding<1>(BASIS_STENCIL_UNIFORM<VECTOR<float,1>,1> const&) const;
template struct BASIS_STENCIL_UNIFORM<VECTOR<float,2>,0>;
template struct BASIS_STENCIL_UNIFORM<VECTOR<float,2>,1>;
template int BASIS_STENCIL_UNIFORM<VECTOR<float,2>,0>::Overlap_Padding<0>(BASIS_STENCIL_UNIFORM<VECTOR<float,2>,0> const&) const;
template int BASIS_STENCIL_UNIFORM<VECTOR<float,2>,0>::Overlap_Padding<1>(BASIS_STENCIL_UNIFORM<VECTOR<float,2>,1> const&) const;
template int BASIS_STENCIL_UNIFORM<VECTOR<float,2>,1>::Overlap_Padding<0>(BASIS_STENCIL_UNIFORM<VECTOR<float,2>,0> const&) const;
template int BASIS_STENCIL_UNIFORM<VECTOR<float,2>,1>::Overlap_Padding<1>(BASIS_STENCIL_UNIFORM<VECTOR<float,2>,1> const&) const;
template struct BASIS_STENCIL_UNIFORM<VECTOR<float,3>,0>;
template struct BASIS_STENCIL_UNIFORM<VECTOR<float,3>,1>;
template int BASIS_STENCIL_UNIFORM<VECTOR<float,3>,0>::Overlap_Padding<0>(BASIS_STENCIL_UNIFORM<VECTOR<float,3>,0> const&) const;
template int BASIS_STENCIL_UNIFORM<VECTOR<float,3>,0>::Overlap_Padding<1>(BASIS_STENCIL_UNIFORM<VECTOR<float,3>,1> const&) const;
template int BASIS_STENCIL_UNIFORM<VECTOR<float,3>,1>::Overlap_Padding<0>(BASIS_STENCIL_UNIFORM<VECTOR<float,3>,0> const&) const;
template int BASIS_STENCIL_UNIFORM<VECTOR<float,3>,1>::Overlap_Padding<1>(BASIS_STENCIL_UNIFORM<VECTOR<float,3>,1> const&) const;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template struct BASIS_STENCIL_UNIFORM<VECTOR<double,1>,0>;
template struct BASIS_STENCIL_UNIFORM<VECTOR<double,1>,1>;
template int BASIS_STENCIL_UNIFORM<VECTOR<double,1>,0>::Overlap_Padding<0>(BASIS_STENCIL_UNIFORM<VECTOR<double,1>,0> const&) const;
template int BASIS_STENCIL_UNIFORM<VECTOR<double,1>,0>::Overlap_Padding<1>(BASIS_STENCIL_UNIFORM<VECTOR<double,1>,1> const&) const;
template int BASIS_STENCIL_UNIFORM<VECTOR<double,1>,1>::Overlap_Padding<0>(BASIS_STENCIL_UNIFORM<VECTOR<double,1>,0> const&) const;
template int BASIS_STENCIL_UNIFORM<VECTOR<double,1>,1>::Overlap_Padding<1>(BASIS_STENCIL_UNIFORM<VECTOR<double,1>,1> const&) const;
template struct BASIS_STENCIL_UNIFORM<VECTOR<double,2>,0>;
template struct BASIS_STENCIL_UNIFORM<VECTOR<double,2>,1>;
template int BASIS_STENCIL_UNIFORM<VECTOR<double,2>,0>::Overlap_Padding<0>(BASIS_STENCIL_UNIFORM<VECTOR<double,2>,0> const&) const;
template int BASIS_STENCIL_UNIFORM<VECTOR<double,2>,0>::Overlap_Padding<1>(BASIS_STENCIL_UNIFORM<VECTOR<double,2>,1> const&) const;
template int BASIS_STENCIL_UNIFORM<VECTOR<double,2>,1>::Overlap_Padding<0>(BASIS_STENCIL_UNIFORM<VECTOR<double,2>,0> const&) const;
template int BASIS_STENCIL_UNIFORM<VECTOR<double,2>,1>::Overlap_Padding<1>(BASIS_STENCIL_UNIFORM<VECTOR<double,2>,1> const&) const;
template struct BASIS_STENCIL_UNIFORM<VECTOR<double,3>,0>;
template struct BASIS_STENCIL_UNIFORM<VECTOR<double,3>,1>;
template int BASIS_STENCIL_UNIFORM<VECTOR<double,3>,0>::Overlap_Padding<0>(BASIS_STENCIL_UNIFORM<VECTOR<double,3>,0> const&) const;
template int BASIS_STENCIL_UNIFORM<VECTOR<double,3>,0>::Overlap_Padding<1>(BASIS_STENCIL_UNIFORM<VECTOR<double,3>,1> const&) const;
template int BASIS_STENCIL_UNIFORM<VECTOR<double,3>,1>::Overlap_Padding<0>(BASIS_STENCIL_UNIFORM<VECTOR<double,3>,0> const&) const;
template int BASIS_STENCIL_UNIFORM<VECTOR<double,3>,1>::Overlap_Padding<1>(BASIS_STENCIL_UNIFORM<VECTOR<double,3>,1> const&) const;
#endif
