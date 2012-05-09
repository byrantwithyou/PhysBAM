//#####################################################################
// Copyright 2012, Craig Schroeder, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Geometry/Basic_Geometry/LINE_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry/PLANE.h>
#include <PhysBAM_Geometry/Basic_Geometry/SEGMENT_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <PhysBAM_Geometry/Finite_Elements/BASIS_INTEGRATION_UNIFORM_NEW.h>
#include <PhysBAM_Geometry/Finite_Elements/CELL_DOMAIN_INTERFACE_NEW.h>
#include <PhysBAM_Geometry/Grids_Uniform_Computations/MARCHING_CUBES.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV,int static_degree> BASIS_INTEGRATION_UNIFORM_NEW<TV,static_degree>::
BASIS_INTEGRATION_UNIFORM_NEW(const GRID<TV>& grid_input,const GRID<TV>& phi_grid_input,const ARRAY<T,TV_INT>& phi_input,CELL_DOMAIN_INTERFACE_NEW<TV>& cdi_input)
    :grid(grid_input),phi_grid(phi_grid_input),cdi(cdi_input),phi(phi_input)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV,int static_degree> BASIS_INTEGRATION_UNIFORM_NEW<TV,static_degree>::
~BASIS_INTEGRATION_UNIFORM_NEW()
{
    volume_blocks.Delete_Pointers_And_Clean_Memory();
    interface_blocks.Delete_Pointers_And_Clean_Memory();
}
//#####################################################################
// Function Precomputed_Integral
//#####################################################################
template<class T,int rank,int sdp1,int d> static T
Precomputed_Integral(const STATIC_TENSOR<T,rank,sdp1>& precompute,const STATIC_POLYNOMIAL<T,rank,d>& poly)
{
    T total=0;
    RANGE<VECTOR<int,rank> > range(VECTOR<int,rank>(),poly.size+1);
    for(RANGE_ITERATOR<rank> it(range);it.Valid();it.Next())
        if(T coeff=poly.terms(it.index))
            total+=coeff*precompute(it.index);
    return total;
}
//#####################################################################
// Function Compute_Entries
//#####################################################################
template<class TV,int static_degree> void BASIS_INTEGRATION_UNIFORM_NEW<TV,static_degree>::
Compute_Entries()
{
    VECTOR<ARRAY<T_FACE>,(1<<TV::m)> sides,interface;
    VECTOR<bool,(1<<TV::m)> enclose_inside;
    VECTOR<int,(1<<TV::m)> direction;
    VECTOR<TV,TV::m> orientation;

    Compute_Open_Entries();

    const VECTOR<TV_INT,(1<<TV::m)>& phi_offsets=GRID<TV>::Binary_Counts(TV_INT());
    const VECTOR<TV_INT,(1<<TV::m)>& phi_double_offsets=GRID<TV>::Binary_Counts(TV_INT(),2);
    const int all_inside=(1<<(1<<TV::m))-1;
    const int all_outside=0;

    int cut_cell_index=-1;

    for(UNIFORM_GRID_ITERATOR_CELL<TV> it(grid);it.Valid();it.Next()){
        int cell_corners=0;
        TV_INT phi_base=it.index*2;
        for(int b=0;b<(1<<TV::m);b++) cell_corners|=(phi(phi_base+phi_double_offsets(b))<0)<<b;
        if(cell_corners==all_inside){Add_Uncut_Cell(it.index,false);continue;}
        if(cell_corners==all_outside){Add_Uncut_Cell(it.index,true);continue;}

        cut_cell_index++;

        for(int b=0;b<(1<<TV::m);b++){
            sides(b).Remove_All();
            interface(b).Remove_All();
            MARCHING_CUBES<TV>::Get_Elements_For_Cell(interface(b),sides(b),direction(b),enclose_inside(b),phi,phi_base+phi_offsets(b));}

        // TODO: Compute weighted normal and set tangent directions here

        for(int b=0;b<(1<<TV::m);b++)
            if(!interface(b).m) Add_Uncut_Fine_Cell(it.index,b,!enclose_inside(b));
            else Add_Cut_Fine_Cell(it.index,b,interface(b),sides(b),direction(b),enclose_inside(b),cut_cell_index,orientation);}
}
//#####################################################################
// Function Compute_Open_Entries
//#####################################################################
template<class TV,int static_degree> void BASIS_INTEGRATION_UNIFORM_NEW<TV,static_degree>::
Compute_Open_Entries()
{
    STATIC_TENSOR<T,TV::m,static_degree+1> uncut_subcell[1<<TV::m];
    const VECTOR<TV_INT,(1<<TV::m)>& counts=GRID<TV>::Binary_Counts(TV_INT());
    T tol=0;

    RANGE<TV_INT> range(TV_INT(),TV_INT()+static_degree+1);
    for(int i=0;i<(1<<TV::m);i++){
        RANGE<TV> subcell_range;
        subcell_range.max_corner=TV(counts(i))*grid.dX/2;
        subcell_range.min_corner=subcell_range.max_corner-grid.dX/2;

        for(RANGE_ITERATOR<TV::m> it(range);it.Valid();it.Next())
            if(volume_monomials_needed(it.index)){
                STATIC_POLYNOMIAL<T,TV::m,static_degree> monomial;
                monomial.Set_Term(it.index,1);
                uncut_subcell[i](it.index)=monomial.Definite_Integral(subcell_range);}
        tol=max(tol,uncut_subcell[i](TV_INT()));}
    tol*=1e-14;

    for(int k=0;k<volume_blocks.m;k++){
        VOLUME_BLOCK* vb=volume_blocks(k);
        for(int i=0;i<vb->overlap_polynomials.m;i++){
            typename VOLUME_BLOCK::OVERLAP_POLYNOMIAL& op=vb->overlap_polynomials(i);
            for(int b=0;b<(1<<TV::m);b++)
                if(op.subcell&(1<<b)){
                    T integral=Precomputed_Integral(uncut_subcell[b],op.polynomial);
                    if(fabs(integral)<tol) continue;
                    OPEN_ENTRY me={op.flat_index_offset,op.flat_index_diff_ref,integral};
                    vb->open_subcell_entries[b].Append(me);}}

        for(int b=0;b<(1<<TV::m);b++){
            vb->open_subcell_entries[b].Coalesce();
            vb->open_entries.Append_Elements(vb->open_subcell_entries[b]);}
        vb->open_entries.Coalesce();}
}
//#####################################################################
// Function Add_Uncut_Cell
//#####################################################################
template<class TV,int static_degree> void BASIS_INTEGRATION_UNIFORM_NEW<TV,static_degree>::
Add_Uncut_Cell(const TV_INT& cell,int enclose_inside)
{
    for(int i=0;i<volume_blocks.m;i++)
        volume_blocks(i)->Add_Open_Entries(cdi.Flatten(cell),enclose_inside);
}
//#####################################################################
// Function Add_Uncut_Fine_Cell
//#####################################################################
template<class TV,int static_degree> void BASIS_INTEGRATION_UNIFORM_NEW<TV,static_degree>::
Add_Uncut_Fine_Cell(const TV_INT& cell,int block,int enclose_inside)
{
    for(int i=0;i<volume_blocks.m;i++)
        volume_blocks(i)->Add_Open_Subcell_Entries(cdi.Flatten(cell),block,enclose_inside);
}
//#####################################################################
// Function Volume
//#####################################################################
template<class T> static T
Volume(VECTOR<VECTOR<T,3>,3> X,int direction)
{
    int sign=direction?-1:1;
    exchange(X(0)(direction),X(0).x);
    exchange(X(1)(direction),X(1).x);
    exchange(X(2)(direction),X(2).x);
    return (X(1).y*(X(2).z-X(0).z)+X(0).y*(X(1).z-X(2).z)+X(2).y*(X(0).z-X(1).z))*(X(0).x+X(1).x+X(2).x)/6*sign;
}
//#####################################################################
// Function Volume
//#####################################################################
template<class T> static T
Volume(VECTOR<VECTOR<T,2>,2> X,int direction)
{
    int sign=direction?-1:1;
    exchange(X(0)(direction),X(0).x);
    exchange(X(1)(direction),X(1).x);
    return (X(1).y-X(0).y)*(X(0).x+X(1).x)/2*sign;
}
//#####################################################################
// Function Add_Cut_Fine_Cell
//#####################################################################
template<class TV,int static_degree> void BASIS_INTEGRATION_UNIFORM_NEW<TV,static_degree>::
Add_Cut_Fine_Cell(const TV_INT& cell,int block,ARRAY<T_FACE>& interface,ARRAY<T_FACE>& sides,
    int direction,bool enclose_inside,int cut_cell_index,VECTOR<TV,TV::m>& orientation)
{
    assert(sides.m);
    assert(interface.m);

    for(int i=0;i<interface.m;i++)
        for(int j=0;j<TV::m;j++)
            interface(i).X(j)=(interface(i).X(j)-(T).5)*grid.dX;

    for(int i=0;i<sides.m;i++)
        for(int j=0;j<TV::m;j++)
            sides(i).X(j)=(sides(i).X(j)-(T).5)*grid.dX;
            
    sides.Append_Elements(interface);

    STATIC_TENSOR<T,TV::m,static_degree+1> precomputed_integrals;
    RANGE<TV_INT> range(TV_INT(),TV_INT()+static_degree+1);
    for(RANGE_ITERATOR<TV::m> it(range);it.Valid();it.Next())
        if(volume_monomials_needed(it.index)){
            STATIC_POLYNOMIAL<T,TV::m,static_degree+1> monomial;
            monomial.Set_Term(it.index,1);
            monomial=monomial.Integrate(direction);
            T integral=0;
            for(int i=0;i<sides.m;i++){
                const VECTOR<TV,TV::m>& V=sides(i).X;
                integral+=monomial.Quadrature_Over_Primitive(V)*T_FACE::Normal(V)(direction);}
            precomputed_integrals(it.index)+=integral;}

    for(int i=0;i<volume_blocks.m;i++){
        VOLUME_BLOCK* vb=volume_blocks(i);
        for(int j=0;j<vb->overlap_polynomials.m;j++){
            typename VOLUME_BLOCK::OVERLAP_POLYNOMIAL& op=vb->overlap_polynomials(j);
            if(op.subcell&(1<<block)){
                T integral=Precomputed_Integral(precomputed_integrals,op.polynomial);
                int flat_index=cdi.Flatten(cell)+op.flat_index_offset;
                vb->Add_Entry(flat_index,op.flat_index_diff_ref,enclose_inside,integral);
                vb->Add_Entry(flat_index,op.flat_index_diff_ref,!enclose_inside,-integral);}}}

    Add_Uncut_Fine_Cell(cell,block,!enclose_inside);

    STATIC_TENSOR<T,TV::m,static_degree+1> precomputed_interface_integrals;
    for(RANGE_ITERATOR<TV::m> it(range);it.Valid();it.Next())
        if(surface_monomials_needed(it.index)){
            STATIC_POLYNOMIAL<T,TV::m,static_degree> monomial;
            monomial.Set_Term(it.index,1);
            for(int i=0;i<interface.m;i++){
                const VECTOR<TV,TV::m> V=interface(i).X;
                precomputed_interface_integrals(it.index)+=monomial.Quadrature_Over_Primitive(V);}}
    
    int sign=enclose_inside?1:-1;

    for(int i=0;i<interface_blocks.m;i++){
        INTERFACE_BLOCK* ib=interface_blocks(i);
        for(int j=0;j<ib->overlap_polynomials.m;j++){
            typename INTERFACE_BLOCK::OVERLAP_POLYNOMIAL& op=ib->overlap_polynomials(j);
            if(op.subcell&(1<<block)){
                T integral=Precomputed_Integral(precomputed_interface_integrals,op.polynomial)*sign;
                ib->Add_Entry(cut_cell_index,op.flat_index_diff_ref,enclose_inside,integral);
                ib->Add_Entry(cut_cell_index,op.flat_index_diff_ref,!enclose_inside,-integral);}}}
}
//#####################################################################
// Function Add_Volume_Block
//#####################################################################
template<class TV,int static_degree> template<int d0,int d1> void BASIS_INTEGRATION_UNIFORM_NEW<TV,static_degree>::
Add_Volume_Block(SYSTEM_VOLUME_BLOCK_HELPER_NEW<TV>& helper,const BASIS_STENCIL_UNIFORM<TV,d0>& s0,
    const BASIS_STENCIL_UNIFORM<TV,d1>& s1,const VECTOR<T,2>& scale)
{
    VOLUME_BLOCK* vb=new VOLUME_BLOCK;
    vb->Initialize(helper,s0,s1,scale);
    volume_blocks.Append(vb);

    for(int i=0;i<vb->overlap_polynomials.m;i++){
        RANGE<TV_INT> range(TV_INT(),vb->overlap_polynomials(i).polynomial.size+1);
        for(RANGE_ITERATOR<TV::m> it(range);it.Valid();it.Next())
            if(vb->overlap_polynomials(i).polynomial.terms(it.index))
                volume_monomials_needed(it.index)=true;}
}
//#####################################################################
// Function Add_Interface_Block
//#####################################################################
template<class TV,int static_degree> template<int d> void BASIS_INTEGRATION_UNIFORM_NEW<TV,static_degree>::
Add_Interface_Block(SYSTEM_INTERFACE_BLOCK_HELPER_NEW<TV>& helper,const BASIS_STENCIL_UNIFORM<TV,d>& s,
    T scale,bool ignore_orientation)
{
    INTERFACE_BLOCK* ib=new INTERFACE_BLOCK;
    ib->Initialize(helper,s,scale,ignore_orientation);
    interface_blocks.Append(ib);
        
    for(int i=0;i<ib->overlap_polynomials.m;i++){
        RANGE<TV_INT> range(TV_INT(),ib->overlap_polynomials(i).polynomial.size+1);
        for(RANGE_ITERATOR<TV::m> it(range);it.Valid();it.Next())
            if(ib->overlap_polynomials(i).polynomial.terms(it.index))
                surface_monomials_needed(it.index)=true;}
}
template class BASIS_INTEGRATION_UNIFORM_NEW<VECTOR<float,3>,2>;
template class BASIS_INTEGRATION_UNIFORM_NEW<VECTOR<float,2>,2>;
template void BASIS_INTEGRATION_UNIFORM_NEW<VECTOR<float,3>,2>::Add_Volume_Block<0,1>(SYSTEM_VOLUME_BLOCK_HELPER_NEW<VECTOR<float,3> >&,
    BASIS_STENCIL_UNIFORM<VECTOR<float,3>,0> const&,BASIS_STENCIL_UNIFORM<VECTOR<float,3>,1> const&,VECTOR<float,2> const&);
template void BASIS_INTEGRATION_UNIFORM_NEW<VECTOR<float,3>,2>::Add_Volume_Block<1,1>(SYSTEM_VOLUME_BLOCK_HELPER_NEW<VECTOR<float,3> >&,
    BASIS_STENCIL_UNIFORM<VECTOR<float,3>,1> const&,BASIS_STENCIL_UNIFORM<VECTOR<float,3>,1> const&,VECTOR<float,2> const&);
template void BASIS_INTEGRATION_UNIFORM_NEW<VECTOR<float,3>,2>::Add_Interface_Block<1>(
    SYSTEM_INTERFACE_BLOCK_HELPER_NEW<VECTOR<float,3> >&,BASIS_STENCIL_UNIFORM<VECTOR<float,3>,1> const&,float,bool);
template void BASIS_INTEGRATION_UNIFORM_NEW<VECTOR<float,2>,2>::Add_Volume_Block<0,1>(SYSTEM_VOLUME_BLOCK_HELPER_NEW<VECTOR<float,2> >&,
    BASIS_STENCIL_UNIFORM<VECTOR<float,2>,0> const&,BASIS_STENCIL_UNIFORM<VECTOR<float,2>,1> const&,VECTOR<float,2> const&);
template void BASIS_INTEGRATION_UNIFORM_NEW<VECTOR<float,2>,2>::Add_Volume_Block<1,1>(SYSTEM_VOLUME_BLOCK_HELPER_NEW<VECTOR<float,2> >&,
    BASIS_STENCIL_UNIFORM<VECTOR<float,2>,1> const&,BASIS_STENCIL_UNIFORM<VECTOR<float,2>,1> const&,VECTOR<float,2> const&);
template void BASIS_INTEGRATION_UNIFORM_NEW<VECTOR<float,2>,2>::Add_Interface_Block<1>(
    SYSTEM_INTERFACE_BLOCK_HELPER_NEW<VECTOR<float,2> >&,BASIS_STENCIL_UNIFORM<VECTOR<float,2>,1> const&,float,bool);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class BASIS_INTEGRATION_UNIFORM_NEW<VECTOR<double,3>,2>;
template class BASIS_INTEGRATION_UNIFORM_NEW<VECTOR<double,2>,2>;
template void BASIS_INTEGRATION_UNIFORM_NEW<VECTOR<double,3>,2>::Add_Volume_Block<0,1>(SYSTEM_VOLUME_BLOCK_HELPER_NEW<VECTOR<double,3> >&,
    BASIS_STENCIL_UNIFORM<VECTOR<double,3>,0> const&,BASIS_STENCIL_UNIFORM<VECTOR<double,3>,1> const&,VECTOR<double,2> const&);
template void BASIS_INTEGRATION_UNIFORM_NEW<VECTOR<double,3>,2>::Add_Volume_Block<1,1>(SYSTEM_VOLUME_BLOCK_HELPER_NEW<VECTOR<double,3> >&,
    BASIS_STENCIL_UNIFORM<VECTOR<double,3>,1> const&,BASIS_STENCIL_UNIFORM<VECTOR<double,3>,1> const&,VECTOR<double,2> const&);
template void BASIS_INTEGRATION_UNIFORM_NEW<VECTOR<double,3>,2>::Add_Interface_Block<1>(
    SYSTEM_INTERFACE_BLOCK_HELPER_NEW<VECTOR<double,3> >&,BASIS_STENCIL_UNIFORM<VECTOR<double,3>,1> const&,double,bool);
template void BASIS_INTEGRATION_UNIFORM_NEW<VECTOR<double,2>,2>::Add_Volume_Block<0,1>(SYSTEM_VOLUME_BLOCK_HELPER_NEW<VECTOR<double,2> >&,
    BASIS_STENCIL_UNIFORM<VECTOR<double,2>,0> const&,BASIS_STENCIL_UNIFORM<VECTOR<double,2>,1> const&,VECTOR<double,2> const&);
template void BASIS_INTEGRATION_UNIFORM_NEW<VECTOR<double,2>,2>::Add_Volume_Block<1,1>(SYSTEM_VOLUME_BLOCK_HELPER_NEW<VECTOR<double,2> >&,
    BASIS_STENCIL_UNIFORM<VECTOR<double,2>,1> const&,BASIS_STENCIL_UNIFORM<VECTOR<double,2>,1> const&,VECTOR<double,2> const&);
template void BASIS_INTEGRATION_UNIFORM_NEW<VECTOR<double,2>,2>::Add_Interface_Block<1>(
    SYSTEM_INTERFACE_BLOCK_HELPER_NEW<VECTOR<double,2> >&,BASIS_STENCIL_UNIFORM<VECTOR<double,2>,1> const&,double,bool);
#endif
