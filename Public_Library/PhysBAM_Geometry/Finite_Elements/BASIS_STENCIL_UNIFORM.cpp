//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform_Arrays/UNIFORM_ARRAY_ITERATOR.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Geometry/Finite_Elements/BASIS_STENCIL_UNIFORM.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> BASIS_STENCIL_UNIFORM<TV>::
BASIS_STENCIL_UNIFORM(TV dx)
    :dX(dx)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> BASIS_STENCIL_UNIFORM<TV>::
~BASIS_STENCIL_UNIFORM()
{
}
//#####################################################################
// Function Add_Symmetric_Entry
//#####################################################################
template<class TV> void BASIS_STENCIL_UNIFORM<TV>::
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
template<class TV> void BASIS_STENCIL_UNIFORM<TV>::
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
template<class TV> void BASIS_STENCIL_UNIFORM<TV>::
Scale_Axes(TV scale)
{
    for(int i=0;i<stencils.m;i++)
        stencils(i).polynomial.Scale_Variables(scale);
}
//#####################################################################
// Function Dice_Stencil
//#####################################################################
template<class TV> void BASIS_STENCIL_UNIFORM<TV>::
Dice_Stencil()
{
    diced.Remove_All();
    const int big_shift=1<<30;
    const int half_big_shift=1<<29;
    RANGE<TV_INT> cell_box=RANGE<TV_INT>::Centered_Box();
    for(int i=0;i<stencils.m;i++){
        RANGE<TV_INT> offset_range=stencils(i).region.Translated(center_offset);
        RANGE<TV_INT> ra((big_shift+1-offset_range.max_corner)/2-half_big_shift,(big_shift+2-offset_range.min_corner)/2-half_big_shift);
        for(UNIFORM_ARRAY_ITERATOR<TV::m> it(ra);it.Valid();it.Next()){
            RANGE<TV_INT> cut_range=RANGE<TV_INT>::Intersect(offset_range.Translated(2*it.index),cell_box);
            DICED e={it.index,cut_range};
            e.polynomial=stencils(i).polynomial;
            e.polynomial.Shift(-dX*(TV(it.index)+TV(center_offset)/2));
            diced.Append(e);}}
}
//#####################################################################
// Function Set_Constant_Stencil
//#####################################################################
template<class TV> void BASIS_STENCIL_UNIFORM<TV>::
Set_Constant_Stencil()
{
    ENTRY e;
    e.region=RANGE<TV_INT>::Centered_Box();
    e.polynomial.terms.Append(MULTIVARIATE_MONOMIAL<TV>(TV_INT(),1));
    e.polynomial.Scale_Variables((T)1/dX);
    stencils.Append(e);
}
//#####################################################################
// Function Set_Multilinear_Stencil
//#####################################################################
template<class TV> void BASIS_STENCIL_UNIFORM<TV>::
Set_Multilinear_Stencil()
{
    ENTRY e;
    e.region.min_corner=TV_INT();
    e.region.max_corner=TV_INT()+2;
    e.polynomial.terms.Append(MULTIVARIATE_MONOMIAL<TV>(TV_INT()+1,TV::m%2==0?1:-1));
    e.polynomial.Shift(TV()-1);
    e.polynomial.Scale_Variables((T)1/dX);
    Add_Symmetric_Entry(e);
}
//#####################################################################
// Function Differentiate
//#####################################################################
template<class TV> void BASIS_STENCIL_UNIFORM<TV>::
Differentiate(int v)
{
    int k=0;
    for(int i=0;i<stencils.m;i++){
        stencils(i).polynomial.Differentiate(v);
        if(stencils(i).polynomial.terms.m)
            stencils(k++)=stencils(i);}
    stencils.Resize(k);
}
//#####################################################################
// Function Print
//#####################################################################
template<class TV> void BASIS_STENCIL_UNIFORM<TV>::
Print() const
{
    LOG::cout<<"Stencils: "<<center_offset<<std::endl;
    for(int i=0;i<stencils.m;i++)
        LOG::cout<<"  "<<stencils(i).region<<"  "<<stencils(i).polynomial<<std::endl;

    LOG::cout<<"Diced:"<<std::endl;
    for(int i=0;i<diced.m;i++)
        LOG::cout<<"  "<<diced(i).index_offset<<"  "<<diced(i).range<<"  "<<diced(i).polynomial<<std::endl;
}
template class BASIS_STENCIL_UNIFORM<VECTOR<float,1> >;
template class BASIS_STENCIL_UNIFORM<VECTOR<float,2> >;
template class BASIS_STENCIL_UNIFORM<VECTOR<float,3> >;
#ifndef COMPILATE_WITHOUT_DOUBLE_SUPPORT
template class BASIS_STENCIL_UNIFORM<VECTOR<double,1> >;
template class BASIS_STENCIL_UNIFORM<VECTOR<double,2> >;
template class BASIS_STENCIL_UNIFORM<VECTOR<double,3> >;
#endif
