//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/UNIFORM_ARRAY_ITERATOR.h>
#include <PhysBAM_Tools/Symbolics/MULTIVARIATE_POLYNOMIAL.h>
#include <PhysBAM_Geometry/Basic_Geometry/LINE_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry/PLANE.h>
#include <PhysBAM_Geometry/Basic_Geometry/SEGMENT_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <PhysBAM_Geometry/Finite_Elements/BASIS_INTEGRATION_UNIFORM.h>
#include <PhysBAM_Geometry/Finite_Elements/BASIS_STENCIL_UNIFORM.h>
#include <PhysBAM_Geometry/Finite_Elements/CELL_MAPPING.h>
#include <PhysBAM_Geometry/Grids_Uniform_Computations/MARCHING_CUBES.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/SYSTEM_MATRIX_HELPER.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> BASIS_INTEGRATION_UNIFORM<TV>::
BASIS_INTEGRATION_UNIFORM(const RANGE<TV_INT> boundary_conditions_input,const GRID<TV>& grid_input,const GRID<TV>& phi_grid_input,
    const BASIS_STENCIL_UNIFORM<TV>& s0_input,const BASIS_STENCIL_UNIFORM<TV>& s1_input,CELL_MAPPING<TV>& cm0_input,
    CELL_MAPPING<TV>& cm1_input,const ARRAY<T,TV_INT>& phi_input)
    :grid(grid_input),phi_grid(phi_grid_input),boundary_conditions(boundary_conditions_input),s0(s0_input),s1(s1_input),
    cm0(cm0_input),cm1(cm1_input),phi(phi_input)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> BASIS_INTEGRATION_UNIFORM<TV>::
~BASIS_INTEGRATION_UNIFORM()
{
}
//#####################################################################
// Function Compute_Overlap_Polynomials
//#####################################################################
template<class TV> void BASIS_INTEGRATION_UNIFORM<TV>::
Compute_Overlap_Polynomials()
{
    for(int i=0;i<s0.diced.m;i++)
        for(int j=0;j<s1.diced.m;j++){
            RANGE<TV_INT> overlap=RANGE<TV_INT>::Intersect(s0.diced(i).range,s1.diced(j).range);
            if(!overlap.Empty_Half_Open()){
                OVERLAP_POLYNOMIALS op = {s0.diced(i).index_offset, s1.diced(j).index_offset, overlap};
                op.polynomial=s0.diced(i).polynomial;
                op.polynomial*=s1.diced(j).polynomial;
                overlap_polynomials.Append(op);}}
}
//#####################################################################
// Function Compute_Open_Entries
//#####################################################################
template<class TV> void BASIS_INTEGRATION_UNIFORM<TV>::
Compute_Open_Entries()
{
    for(int i=0;i<overlap_polynomials.m;i++){
        T integral=overlap_polynomials(i).polynomial.Definite_Integral(RANGE<TV>(overlap_polynomials(i).range)*(grid.dX/2));
        MATRIX_ENTRY me={overlap_polynomials(i).index_offset0,overlap_polynomials(i).index_offset1,integral};
        raw_open_entries.Append(me);}
    open_entries=raw_open_entries;
    open_entries.Coalesce();
}
//#####################################################################
// Function Compute_Matrix
//#####################################################################
template<class TV> void BASIS_INTEGRATION_UNIFORM<TV>::
Compute_Matrix(SYSTEM_MATRIX_HELPER<T>& helper)
{
    Compute_Overlap_Polynomials();
    Compute_Open_Entries();

    int coarse_factor=grid.counts.x/phi_grid.counts.x;
    RANGE<TV_INT> double_coarse_range(TV_INT(),TV_INT()+2*coarse_factor),coarse_range(TV_INT(),TV_INT()+coarse_factor);
    ARRAY<ARRAY<T_FACE>,TV_INT> cut_elements(double_coarse_range);

    helper.New_Block();
    for(UNIFORM_GRID_ITERATOR_CELL<TV> it(phi_grid);it.Valid();it.Next()){
        ARRAY<T_FACE> elements;
        int dir=0;
        bool enclose_inside=true;
        MARCHING_CUBES<TV>::Get_Elements_For_Cell(elements,elements,dir,enclose_inside,phi,it.index);

        if(!elements.m){ // Uncut cell; emit the standard stencil
            for(UNIFORM_ARRAY_ITERATOR<TV::m> it2(coarse_range);it2.Valid();it2.Next())
                Add_Uncut_Stencil(helper,it.index*coarse_factor+it2.index,!enclose_inside);
            continue;}

        RANGE<TV_INT> flat_range(double_coarse_range);
        flat_range.max_corner(dir)=1;
        Cut_Elements(cut_elements,elements,double_coarse_range,RANGE<TV>::Unit_Box(),dir);

        for(UNIFORM_ARRAY_ITERATOR<TV::m> it2(coarse_range);it2.Valid();it2.Next()){
            for(int i=0;i<overlap_polynomials.m;i++){
                ARRAY<T_FACE> elements;
                RANGE<TV_INT> range=overlap_polynomials(i).range;
                range.min_corner(dir)=0;
                range.max_corner(dir)=1;
                for(UNIFORM_ARRAY_ITERATOR<TV::m> it3(range);it3.Valid();it3.Next()){
                    TV_INT index=it2.index*coarse_factor+it3.index+1;
                    index(dir)=0;
                    elements.Append_Elements(cut_elements(index));}
                Add_Cut_Stencil(helper,elements,it.index*coarse_factor+it2.index,dir,enclose_inside,it2.index,i);}}
        for(UNIFORM_ARRAY_ITERATOR<TV::m> it2(flat_range);it2.Valid();it2.Next())
            cut_elements(it2.index).Remove_All();}
}
//#####################################################################
// Function Cut_Elements
//#####################################################################
template<class TV> void BASIS_INTEGRATION_UNIFORM<TV>::
Apply_Matrix_Entry(SYSTEM_MATRIX_HELPER<T>& helper,const TV_INT& cell,bool inside,MATRIX_ENTRY me)
{
    me.index0+=cell;
    me.index1+=cell;
    int i0=cm0.Get_Index(me.index0,inside);
    int i1=cm1.Get_Index(me.index1,inside);
    assert(i0>=0 && i1>=0);
    helper.data.Append(TRIPLE<int,int,T>(i0,i1,me.x));
}
//#####################################################################
// Function Cut_Elements
//#####################################################################
template<class TV> void BASIS_INTEGRATION_UNIFORM<TV>::
Add_Uncut_Stencil(SYSTEM_MATRIX_HELPER<T>& helper,const TV_INT& cell,bool inside)
{
    for(int i=0;i<open_entries.m;i++) Apply_Matrix_Entry(helper,cell,inside,open_entries(i));
}
template<class T> static T
Volume(VECTOR<VECTOR<T,3>,3> X,int dir)
{
    int sign=dir?-1:1;
    exchange(X(0)(dir),X(0).x);
    exchange(X(1)(dir),X(1).x);
    exchange(X(2)(dir),X(2).x);
    return (X(1).y*(X(2).z-X(0).z)+X(0).y*(X(1).z-X(2).z)+X(2).y*(X(0).z-X(1).z))*(X(0).x+X(1).x+X(2).x)/6*sign;
}
template<class T> static T
Volume(VECTOR<VECTOR<T,2>,2> X,int dir)
{
    int sign=dir?-1:1;
    exchange(X(0)(dir),X(0).x);
    exchange(X(1)(dir),X(1).x);
    return (X(1).y-X(0).y)*(X(0).x+X(1).x)/2*sign;
}
//#####################################################################
// Function Cut_Elements
//#####################################################################
template<class TV> void BASIS_INTEGRATION_UNIFORM<TV>::
Add_Cut_Stencil(SYSTEM_MATRIX_HELPER<T>& helper,const ARRAY<T_FACE>& elements,const TV_INT& cell,int dir,bool enclose_inside,const TV_INT& sub_cell,int opi)
{
    if(!elements.m) return Apply_Matrix_Entry(helper,cell,!enclose_inside,raw_open_entries(opi));
    int coarse_factor=grid.counts.x/phi_grid.counts.x;

    T volume_inside=0;
    ARRAY<T_FACE> projected_elements(elements);
    const OVERLAP_POLYNOMIALS& op=overlap_polynomials(opi);
    T mn=(T)op.range.min_corner(dir)/2,mx=(T)op.range.max_corner(dir)/2;
    for(int i=0;i<projected_elements.m;i++){
        for(int j=0;j<TV::m;j++){
            projected_elements(i).X(j)=projected_elements(i).X(j)*coarse_factor-TV(sub_cell)-(T).5;
            projected_elements(i).X(j)(dir)=clamp(projected_elements(i).X(j)(dir),mn,mx);}
        volume_inside+=Volume(reinterpret_cast<const VECTOR<TV,TV::m>&>(projected_elements(i).X(0)),dir);}

    // Check for full or empty cell.
    T full_volume=(T)op.range.Size()/(1<<TV::m);
    T frac=volume_inside/full_volume;

    if(1-frac<1e-14) return Apply_Matrix_Entry(helper,cell,enclose_inside,raw_open_entries(opi));
    Apply_Matrix_Entry(helper,cell,!enclose_inside,raw_open_entries(opi));
    if(frac<1e-14) return;

    MULTIVARIATE_POLYNOMIAL<TV> poly(op.polynomial);
    poly.Integrate(dir);
    TV_INT index0=op.index_offset0+cell;
    TV_INT index1=op.index_offset1+cell;

    T integral=0;
    for(int i=0;i<projected_elements.m;i++){
        for(int j=0;j<TV::m;j++) projected_elements(i).X(j)*=grid.dX;
        const VECTOR<TV,TV::m>& V=reinterpret_cast<const VECTOR<TV,TV::m>&>(projected_elements(i).x1);
        integral+=poly.Integrate_Over_Primitive(V)*T_FACE::Normal(V)(dir);}

    helper.data.Append(TRIPLE<int,int,T>(cm0.Get_Index(index0,enclose_inside),cm1.Get_Index(index1,enclose_inside),integral));
    helper.data.Append(TRIPLE<int,int,T>(cm0.Get_Index(index0,!enclose_inside),cm1.Get_Index(index1,!enclose_inside),-integral));
}
//#####################################################################
// Function Cut_Elements
//#####################################################################
template<class TV> void BASIS_INTEGRATION_UNIFORM<TV>::
Cut_Elements(ARRAY<ARRAY<T_FACE>,TV_INT>& cut_elements,const ARRAY<T_FACE>& elements,const RANGE<TV_INT>& range,const RANGE<TV>& domain,int dir)
{
    TV_INT size=range.Edge_Lengths();
    for(int a=0;a<TV::m;a++)
        if(size(a)>1){
            TV_INT new_size=size;
            new_size(a)/=2;
            ARRAY<T_FACE> t0,t1;
            TV pt=domain.min_corner+domain.Edge_Lengths()*TV(new_size)/TV(size);
            typename BASIC_GEOMETRY_POLICY<TV>::HYPERPLANE plane(TV::Axis_Vector(a),pt);
            for(int i=0;i<elements.m;i++)
                T_FACE::Cut_With_Hyperplane(elements(i),plane,t0,t1,1e-14);
            RANGE<TV> domain0(domain),domain1(domain);
            RANGE<TV_INT> range0(range),range1(range);
            domain0.max_corner(a)=domain1.min_corner(a)=pt(a);
            range0.max_corner(a)=range1.min_corner(a)=range.min_corner(a)+new_size(a);
            Cut_Elements(cut_elements,t0,range0,domain0,dir);
            Cut_Elements(cut_elements,t1,range1,domain1,dir);
            return;}
    TV_INT cell=range.min_corner;
    cell(dir)=0;
    cut_elements(cell).Append_Elements(elements);
}
template class BASIS_INTEGRATION_UNIFORM<VECTOR<float,2> >;
template class BASIS_INTEGRATION_UNIFORM<VECTOR<float,3> >;
#ifndef COMPILATE_WITHOUT_DOUBLE_SUPPORT
template class BASIS_INTEGRATION_UNIFORM<VECTOR<double,2> >;
template class BASIS_INTEGRATION_UNIFORM<VECTOR<double,3> >;
#endif
