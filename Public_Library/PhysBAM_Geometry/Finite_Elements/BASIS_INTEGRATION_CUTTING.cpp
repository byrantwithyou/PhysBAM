//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Symbolics/STATIC_POLYNOMIAL.h>
#include <PhysBAM_Geometry/Basic_Geometry/LINE_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry/PLANE.h>
#include <PhysBAM_Geometry/Basic_Geometry/SEGMENT_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <PhysBAM_Geometry/Finite_Elements/BASIS_INTEGRATION_CUTTING.h>
#include <PhysBAM_Geometry/Finite_Elements/BASIS_STENCIL_UNIFORM.h>
#include <PhysBAM_Geometry/Finite_Elements/CELL_MAPPING.h>
#include <PhysBAM_Geometry/Grids_Uniform_Computations/MARCHING_CUBES.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/SYSTEM_MATRIX_HELPER.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV,int static_degree> BASIS_INTEGRATION_CUTTING<TV,static_degree>::
BASIS_INTEGRATION_CUTTING(const RANGE<TV_INT> boundary_conditions_input,const GRID<TV>& grid_input,const GRID<TV>& phi_grid_input,
    const ARRAY<T,TV_INT>& phi_input)
    :grid(grid_input),phi_grid(phi_grid_input),boundary_conditions(boundary_conditions_input),phi(phi_input),
    coarse_factor(grid.counts.x/phi_grid.counts.x),double_coarse_range(TV_INT(),TV_INT()+2*coarse_factor),
    coarse_range(TV_INT(),TV_INT()+coarse_factor)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV,int static_degree> BASIS_INTEGRATION_CUTTING<TV,static_degree>::
~BASIS_INTEGRATION_CUTTING()
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
// Function Compute
//#####################################################################
template<class TV,int static_degree> void BASIS_INTEGRATION_CUTTING<TV,static_degree>::
Compute()
{
    ARRAY<ARRAY<PAIR<T_FACE,int> >,TV_INT> cut_sides(double_coarse_range),cut_interface(double_coarse_range);
    ARRAY<T_FACE> sides,interface,array_one(1);

    Compute_Open_Entries();
    const VECTOR<TV_INT,(1<<TV::m)>& counts=GRID<TV>::Binary_Counts(TV_INT());

    int element_base=0;
    for(UNIFORM_GRID_ITERATOR_CELL<TV> it(phi_grid);it.Valid();it.Next()){
        sides.Remove_All();
        interface.Remove_All();
        int dir=0;
        bool enclose_inside=true;
        MARCHING_CUBES<TV>::Get_Elements_For_Cell(interface,sides,dir,enclose_inside,phi,it.index);

        if(!interface.m){
            Add_Uncut_Coarse_Cell(it.index,!enclose_inside);
            continue;}

        RANGE<TV_INT> flat_range(double_coarse_range);
        flat_range.max_corner(dir)=1;
        for(int e=0;e<interface.m;e++){
            array_one(0)=interface(e);
            Cut_Elements(cut_interface,array_one,double_coarse_range,RANGE<TV>::Unit_Box(),-1,e);}
        Cut_Elements(cut_sides,sides,double_coarse_range,RANGE<TV>::Unit_Box(),dir,-1);

        for(RANGE_ITERATOR<TV::m> it2(double_coarse_range);it2.Valid();it2.Next()){
            TV_INT side_index=it2.index;
            side_index(dir)=0;
            cut_sides(side_index).Append_Elements(cut_interface(it2.index));}

        for(RANGE_ITERATOR<TV::m> it2(coarse_range);it2.Valid();it2.Next()){
            TV_INT cell_index=it.index*coarse_factor+it2.index;
            for(int b=0;b<(1<<TV::m);b++){
                TV_INT side_index=it2.index*2+counts(b);
                ARRAY<PAIR<T_FACE,int> >& interface_elements=cut_interface(side_index);
                side_index(dir)=0;
                Add_Cut_Subcell(cut_sides(side_index),interface_elements,cell_index,it2.index,dir,enclose_inside,b,element_base);
                interface_elements.Remove_All();}}
        element_base+=interface.m;

        for(RANGE_ITERATOR<TV::m> it2(double_coarse_range);it2.Valid();it2.Next())
            cut_interface(it2.index).Remove_All();
        for(RANGE_ITERATOR<TV::m> it2(flat_range);it2.Valid();it2.Next())
            cut_sides(it2.index).Remove_All();}
}
//#####################################################################
// Function Compute_Open_Entries
//#####################################################################
template<class TV,int static_degree> void BASIS_INTEGRATION_CUTTING<TV,static_degree>::
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
            if(monomials_needed(it.index)){
                STATIC_POLYNOMIAL<T,TV::m,static_degree> monomial;
                monomial.Set_Term(it.index,1);
                uncut_subcell[i](it.index)=monomial.Definite_Integral(subcell_range);}
        tol=max(tol,uncut_subcell[i](TV_INT()));}
    tol*=1e-14;

    for(int k=0;k<volume_blocks.m;k++){
        VOLUME_BLOCK* vb=volume_blocks(k);
        for(int i=0;i<vb->overlap.m;i++){
            for(int b=0;b<(1<<TV::m);b++)
                if(vb->overlap(i).subcell&(1<<b)){
                    T integral=Precomputed_Integral(uncut_subcell[b],vb->overlap(i).polynomial);
                    if(fabs(integral)<tol) continue;
                    VOLUME_MATRIX_ENTRY me={vb->overlap(i).index_offset0,vb->overlap(i).index_offset1,integral};
                    vb->open_subcell_entries[b].Append(me);}}

        for(int b=0;b<(1<<TV::m);b++){
            vb->open_subcell_entries[b].Coalesce();
            vb->open_entries.Append_Elements(vb->open_subcell_entries[b]);}
        vb->open_entries.Coalesce();}
}
//#####################################################################
// Function Apply_Matrix_Entry
//#####################################################################
template<class TV,int static_degree> void BASIS_INTEGRATION_CUTTING<TV,static_degree>::
Apply_Matrix_Entry(VOLUME_BLOCK* vb,const TV_INT& cell,bool inside,VOLUME_MATRIX_ENTRY me)
{
    me.index0+=cell;
    me.index1+=cell;
    int i0=vb->cm0->Get_Index(me.index0,inside);
    int i1=vb->cm1->Get_Index(me.index1,inside);
    assert(i0>=0 && i1>=0);
    vb->helper->data.Append(TRIPLE<int,int,T>(i0,i1,me.x*vb->scale));
}
//#####################################################################
// Function Add_Uncut_Cell
//#####################################################################
template<class TV,int static_degree> void BASIS_INTEGRATION_CUTTING<TV,static_degree>::
Add_Uncut_Coarse_Cell(const TV_INT& coarse_cell,int inside)
{
    for(RANGE_ITERATOR<TV::m> it2(coarse_range);it2.Valid();it2.Next()){
        TV_INT index=coarse_cell*coarse_factor+it2.index;
        for(int i=0;i<volume_blocks.m;i++){
            VOLUME_BLOCK* vb=volume_blocks(i);
            for(int j=0;j<vb->open_entries.m;j++)
                Apply_Matrix_Entry(vb,index,inside,vb->open_entries(j));}}
}
//#####################################################################
// Function Add_Uncut_Cell
//#####################################################################
template<class TV,int static_degree> void BASIS_INTEGRATION_CUTTING<TV,static_degree>::
Add_Uncut_Fine_Cell(const TV_INT& cell,int block,int inside)
{
    for(int i=0;i<volume_blocks.m;i++){
        VOLUME_BLOCK* vb=volume_blocks(i);
        for(int j=0;j<vb->open_subcell_entries[block].m;j++)
            Apply_Matrix_Entry(vb,cell,inside,vb->open_subcell_entries[block](j));}
}
//#####################################################################
// Function Cut_Elements
//#####################################################################
template<class TV,int static_degree> void BASIS_INTEGRATION_CUTTING<TV,static_degree>::
Cut_Elements(ARRAY<ARRAY<PAIR<T_FACE,int> >,TV_INT>& cut_elements,const ARRAY<T_FACE>& elements,const RANGE<TV_INT>& range,
    const RANGE<TV>& domain,int dir,int e)
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
            Cut_Elements(cut_elements,t0,range0,domain0,dir,e);
            Cut_Elements(cut_elements,t1,range1,domain1,dir,e);
            return;}
    TV_INT cell=range.min_corner;
    if(dir>=0) cell(dir)=0;
    
    for(int i=0;i<elements.m;i++)
        cut_elements(cell).Append(PAIR<T_FACE,int>(elements(i),e));
}
//#####################################################################
// Function Add_Block
//#####################################################################
template<class TV,int static_degree> template<int d0,int d1> int BASIS_INTEGRATION_CUTTING<TV,static_degree>::
Add_Block(SYSTEM_MATRIX_HELPER<T>& helper,const BASIS_STENCIL_UNIFORM<TV,d0>& s0,const BASIS_STENCIL_UNIFORM<TV,d1>& s1,
    CELL_MAPPING<TV>& cm0,CELL_MAPPING<TV>& cm1,T scale)
{
    VOLUME_BLOCK* vb=new VOLUME_BLOCK;
    vb->cm0=&cm0;
    vb->cm1=&cm1;
    vb->scale=scale;
    vb->helper=&helper;
    for(int i=0;i<s0.diced.m;i++)
        for(int j=0;j<s1.diced.m;j++){
            int overlap=s0.diced(i).subcell&s1.diced(j).subcell;
            if(overlap){
                OVERLAP_VOLUME_POLYNOMIALS op = {s0.diced(i).index_offset, s1.diced(j).index_offset, overlap};
                op.polynomial=s0.diced(i).polynomial*s1.diced(j).polynomial;
                vb->overlap.Append(op);}}

    for(int i=0;i<vb->overlap.m;i++){
        RANGE<TV_INT> range(TV_INT(),vb->overlap(i).polynomial.size+1);
        for(RANGE_ITERATOR<TV::m> it(range);it.Valid();it.Next())
            if(vb->overlap(i).polynomial.terms(it.index))
                monomials_needed(it.index)=true;}

    return volume_blocks.Append(vb);
}
//#####################################################################
// Function Add_Block
//#####################################################################
template<class TV,int static_degree> template<int d> int BASIS_INTEGRATION_CUTTING<TV,static_degree>::
Add_Block(SYSTEM_MATRIX_HELPER<T>& helper,const BASIS_STENCIL_UNIFORM<TV,d>& s,CELL_MAPPING<TV>& cm,T scale)
{
    INTERFACE_BLOCK* ib=new INTERFACE_BLOCK;
    ib->cm=&cm;
    ib->overlap.Resize(s.diced.m);
    for(int i=0;i<ib->overlap.m;i++){
        ib->overlap(i).index_offset=s.diced(i).index_offset;
        ib->overlap(i).subcell=s.diced(i).subcell;
        ib->overlap(i).polynomial=s.diced(i).polynomial;}
    ib->scale=scale;
    ib->helper=&helper;
    for(int i=0;i<ib->overlap.m;i++){
        RANGE<TV_INT> range(TV_INT(),ib->overlap(i).polynomial.size+1);
        for(RANGE_ITERATOR<TV::m> it(range);it.Valid();it.Next())
            if(ib->overlap(i).polynomial.terms(it.index))
                monomials_needed(it.index)=true;}

    return interface_blocks.Append(ib);
}
//#####################################################################
// Function Volume
//#####################################################################
template<class T> static T
Volume(VECTOR<VECTOR<T,3>,3> X,int dir)
{
    int sign=dir?-1:1;
    exchange(X(0)(dir),X(0).x);
    exchange(X(1)(dir),X(1).x);
    exchange(X(2)(dir),X(2).x);
    return (X(1).y*(X(2).z-X(0).z)+X(0).y*(X(1).z-X(2).z)+X(2).y*(X(0).z-X(1).z))*(X(0).x+X(1).x+X(2).x)/6*sign;
}
//#####################################################################
// Function Volume
//#####################################################################
template<class T> static T
Volume(VECTOR<VECTOR<T,2>,2> X,int dir)
{
    int sign=dir?-1:1;
    exchange(X(0)(dir),X(0).x);
    exchange(X(1)(dir),X(1).x);
    return (X(1).y-X(0).y)*(X(0).x+X(1).x)/2*sign;
}
//#####################################################################
// Function Add_Cut_Stencil
//#####################################################################
template<class TV,int static_degree> void BASIS_INTEGRATION_CUTTING<TV,static_degree>::
Add_Cut_Subcell(const ARRAY<PAIR<T_FACE,int> >& side_elements,const ARRAY<PAIR<T_FACE,int> >& interface_elements,
    const TV_INT& cell,const TV_INT& subcell_cell,int dir,bool enclose_inside,int block,int element_base)
{
    const VECTOR<TV_INT,(1<<TV::m)>& counts=GRID<TV>::Binary_Counts(TV_INT());
    if(!side_elements.m){
        Add_Uncut_Fine_Cell(cell,block,!enclose_inside);
        assert(!interface_elements.m);
        return;}

    RANGE<TV> subcell_range;
    subcell_range.max_corner=TV(counts(block))*(grid.dX/2);
    subcell_range.min_corner=subcell_range.max_corner-(grid.dX/2);

    ARRAY<PAIR<T_FACE,int> > projected_elements(side_elements);
    T mn=subcell_range.min_corner(dir),mx=subcell_range.max_corner(dir);
    for(int i=0;i<projected_elements.m;i++){
        for(int j=0;j<TV::m;j++){
            projected_elements(i).x.X(j)=(projected_elements(i).x.X(j)*coarse_factor-TV(subcell_cell)-(T).5)*grid.dX;
            projected_elements(i).x.X(j)(dir)=clamp(projected_elements(i).x.X(j)(dir),mn,mx);}}

    if(!interface_elements.m){
        T volume_inside=0;
        for(int i=0;i<projected_elements.m;i++)
            volume_inside+=Volume(reinterpret_cast<const VECTOR<TV,TV::m>&>(projected_elements(i).x.X(0)),dir);
        bool filled=volume_inside>subcell_range.Size()/2;
        Add_Uncut_Fine_Cell(cell,block,enclose_inside==filled);
        return;}

    STATIC_TENSOR<T,TV::m,static_degree+1> precomputed_integrals;
    RANGE<TV_INT> range(TV_INT(),TV_INT()+static_degree+1);
    for(RANGE_ITERATOR<TV::m> it(range);it.Valid();it.Next())
        if(monomials_needed(it.index)){
            STATIC_POLYNOMIAL<T,TV::m,static_degree+1> monomial;
            monomial.Set_Term(it.index,1);
            monomial=monomial.Integrate(dir);
            T integral=0;
            for(int i=0;i<projected_elements.m;i++){
                const VECTOR<TV,TV::m>& V=reinterpret_cast<const VECTOR<TV,TV::m>&>(projected_elements(i).x.x1);
                integral+=monomial.Integrate_Over_Primitive(V)*T_FACE::Normal(V)(dir);}
            precomputed_integrals(it.index)=integral;}

    Add_Uncut_Fine_Cell(cell,block,!enclose_inside);

    for(int i=0;i<volume_blocks.m;i++){
        VOLUME_BLOCK* vb=volume_blocks(i);
        for(int j=0;j<vb->overlap.m;j++){
            if(vb->overlap(j).subcell&(1<<block)){
                T integral=Precomputed_Integral(precomputed_integrals,vb->overlap(j).polynomial)*vb->scale;
                TV_INT index0=vb->overlap(j).index_offset0+cell;
                TV_INT index1=vb->overlap(j).index_offset1+cell;
                int index_i0=vb->cm0->Get_Index(index0,enclose_inside);
                int index_o0=vb->cm0->Get_Index(index0,!enclose_inside);
                int index_i1=vb->cm1->Get_Index(index1,enclose_inside);
                int index_o1=vb->cm1->Get_Index(index1,!enclose_inside);
                vb->helper->data.Append(TRIPLE<int,int,T>(index_i0,index_i1,integral));
                vb->helper->data.Append(TRIPLE<int,int,T>(index_o0,index_o1,-integral));}}}

    STATIC_TENSOR<T,TV::m,static_degree+1> precomputed_interface_integrals[subcell_elements];
    bool has_element[subcell_elements]={};
    for(RANGE_ITERATOR<TV::m> it(range);it.Valid();it.Next())
        if(monomials_needed(it.index)){
            STATIC_POLYNOMIAL<T,TV::m,static_degree> monomial;
            monomial.Set_Term(it.index,1);
            for(int i=0;i<interface_elements.m;i++){
                VECTOR<TV,TV::m> V=reinterpret_cast<const VECTOR<TV,TV::m>&>(interface_elements(i).x.x1);
                for(int j=0;j<TV::m;j++) V(j)=(V(j)*coarse_factor-TV(counts(block))-(T).5)*grid.dX;
                int e=interface_elements(i).y;
                has_element[e]=true;
                precomputed_interface_integrals[e](it.index)=monomial.Integrate_Over_Primitive(V);}}

    int sign=enclose_inside?1:-1;
    for(int i=0;i<interface_blocks.m;i++){
        INTERFACE_BLOCK* ib=interface_blocks(i);
        for(int j=0;j<ib->overlap.m;j++){
            if(ib->overlap(j).subcell&(1<<block))
                for(int k=0;k<subcell_elements;k++)
                    if(has_element[k]){
                        T integral=Precomputed_Integral(precomputed_interface_integrals[k],ib->overlap(j).polynomial)*ib->scale*sign;
                        TV_INT index=ib->overlap(j).index_offset+cell;
                        int index_i=ib->cm->Get_Index(index,enclose_inside);
                        int index_o=ib->cm->Get_Index(index,!enclose_inside);
                        ib->helper->data.Append(TRIPLE<int,int,T>(index_i,element_base+k,integral));
                        ib->helper->data.Append(TRIPLE<int,int,T>(index_o,element_base+k,-integral));}}}
}
template class BASIS_INTEGRATION_CUTTING<VECTOR<float,2>,2>;
template class BASIS_INTEGRATION_CUTTING<VECTOR<float,3>,2>;
template int BASIS_INTEGRATION_CUTTING<VECTOR<float,3>,2>::Add_Block<1,0>(SYSTEM_MATRIX_HELPER<float>&,
    BASIS_STENCIL_UNIFORM<VECTOR<float,3>,1> const&,BASIS_STENCIL_UNIFORM<VECTOR<float,3>,0> const&,
    CELL_MAPPING<VECTOR<float,3> >&,CELL_MAPPING<VECTOR<float,3> >&,float);
template int BASIS_INTEGRATION_CUTTING<VECTOR<float,3>,2>::Add_Block<1,1>(SYSTEM_MATRIX_HELPER<float>&,
    BASIS_STENCIL_UNIFORM<VECTOR<float,3>,1> const&,BASIS_STENCIL_UNIFORM<VECTOR<float,3>,1> const&,
    CELL_MAPPING<VECTOR<float,3> >&,CELL_MAPPING<VECTOR<float,3> >&,float);
template int BASIS_INTEGRATION_CUTTING<VECTOR<float,3>,2>::Add_Block<1>(SYSTEM_MATRIX_HELPER<float>&,
    BASIS_STENCIL_UNIFORM<VECTOR<float,3>,1> const&,CELL_MAPPING<VECTOR<float,3> >&,float);
template int BASIS_INTEGRATION_CUTTING<VECTOR<float,2>,2>::Add_Block<1,0>(SYSTEM_MATRIX_HELPER<float>&,
    BASIS_STENCIL_UNIFORM<VECTOR<float,2>,1> const&,BASIS_STENCIL_UNIFORM<VECTOR<float,2>,0> const&,
    CELL_MAPPING<VECTOR<float,2> >&,CELL_MAPPING<VECTOR<float,2> >&,float);
template int BASIS_INTEGRATION_CUTTING<VECTOR<float,2>,2>::Add_Block<1,1>(SYSTEM_MATRIX_HELPER<float>&,
    BASIS_STENCIL_UNIFORM<VECTOR<float,2>,1> const&,BASIS_STENCIL_UNIFORM<VECTOR<float,2>,1> const&,
    CELL_MAPPING<VECTOR<float,2> >&,CELL_MAPPING<VECTOR<float,2> >&,float);
template int BASIS_INTEGRATION_CUTTING<VECTOR<float,2>,2>::Add_Block<1>(SYSTEM_MATRIX_HELPER<float>&,
    BASIS_STENCIL_UNIFORM<VECTOR<float,2>,1> const&,CELL_MAPPING<VECTOR<float,2> >&,float);
#ifndef COMPILATE_WITHOUT_DOUBLE_SUPPORT
template class BASIS_INTEGRATION_CUTTING<VECTOR<double,2>,2>;
template class BASIS_INTEGRATION_CUTTING<VECTOR<double,3>,2>;
template int BASIS_INTEGRATION_CUTTING<VECTOR<double,3>,2>::Add_Block<1,0>(SYSTEM_MATRIX_HELPER<double>&,
    BASIS_STENCIL_UNIFORM<VECTOR<double,3>,1> const&,BASIS_STENCIL_UNIFORM<VECTOR<double,3>,0> const&,
    CELL_MAPPING<VECTOR<double,3> >&,CELL_MAPPING<VECTOR<double,3> >&,double);
template int BASIS_INTEGRATION_CUTTING<VECTOR<double,3>,2>::Add_Block<1,1>(SYSTEM_MATRIX_HELPER<double>&,
    BASIS_STENCIL_UNIFORM<VECTOR<double,3>,1> const&,BASIS_STENCIL_UNIFORM<VECTOR<double,3>,1> const&,
    CELL_MAPPING<VECTOR<double,3> >&,CELL_MAPPING<VECTOR<double,3> >&,double);
template int BASIS_INTEGRATION_CUTTING<VECTOR<double,3>,2>::Add_Block<1>(SYSTEM_MATRIX_HELPER<double>&,
    BASIS_STENCIL_UNIFORM<VECTOR<double,3>,1> const&,CELL_MAPPING<VECTOR<double,3> >&,double);
template int BASIS_INTEGRATION_CUTTING<VECTOR<double,2>,2>::Add_Block<1,0>(SYSTEM_MATRIX_HELPER<double>&,
    BASIS_STENCIL_UNIFORM<VECTOR<double,2>,1> const&,BASIS_STENCIL_UNIFORM<VECTOR<double,2>,0> const&,
    CELL_MAPPING<VECTOR<double,2> >&,CELL_MAPPING<VECTOR<double,2> >&,double);
template int BASIS_INTEGRATION_CUTTING<VECTOR<double,2>,2>::Add_Block<1,1>(SYSTEM_MATRIX_HELPER<double>&,
    BASIS_STENCIL_UNIFORM<VECTOR<double,2>,1> const&,BASIS_STENCIL_UNIFORM<VECTOR<double,2>,1> const&,
    CELL_MAPPING<VECTOR<double,2> >&,CELL_MAPPING<VECTOR<double,2> >&,double);
template int BASIS_INTEGRATION_CUTTING<VECTOR<double,2>,2>::Add_Block<1>(SYSTEM_MATRIX_HELPER<double>&,
    BASIS_STENCIL_UNIFORM<VECTOR<double,2>,1> const&,CELL_MAPPING<VECTOR<double,2> >&,double);
#endif
