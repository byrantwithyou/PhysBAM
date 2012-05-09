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
    :grid(grid_input),phi_grid(phi_grid_input),cdi(cdi_input),phi(phi_input),coarse_factor(grid.counts.x/phi_grid.counts.x),
    double_coarse_range(TV_INT(),TV_INT()+2*coarse_factor),coarse_range(TV_INT(),TV_INT()+coarse_factor)
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

        TV_INT base_index=it.index*coarse_factor;
        for(RANGE_ITERATOR<TV::m> it2(coarse_range);it2.Valid();it2.Next()){
            TV_INT cell_index=base_index+it2.index;
            for(int b=0;b<(1<<TV::m);b++){
                TV_INT side_index=it2.index*2+counts(b);
                ARRAY<PAIR<T_FACE,int> >& interface_elements=cut_interface(side_index);
                side_index(dir)=0;
                Add_Cut_Subcell(cut_sides(side_index),interface_elements,cell_index,it2.index,dir,enclose_inside,b,element_base);
                interface_elements.Remove_All();}}
        cdi.Set_Flat_Base(element_base,base_index);
        element_base+=interface.m;

        for(RANGE_ITERATOR<TV::m> it2(double_coarse_range);it2.Valid();it2.Next())
            cut_interface(it2.index).Remove_All();
        for(RANGE_ITERATOR<TV::m> it2(flat_range);it2.Valid();it2.Next())
            cut_sides(it2.index).Remove_All();}
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
// Function Add_Uncut_Coarse_Cell
//#####################################################################
template<class TV,int static_degree> void BASIS_INTEGRATION_UNIFORM_NEW<TV,static_degree>::
Add_Uncut_Coarse_Cell(const TV_INT& coarse_cell,int inside)
{
    for(RANGE_ITERATOR<TV::m> it2(coarse_range);it2.Valid();it2.Next()){
        TV_INT cell=coarse_cell*coarse_factor+it2.index;
        for(int i=0;i<volume_blocks.m;i++)
            volume_blocks(i)->Add_Open_Entries(cdi.Flatten(cell),inside);}
}
//#####################################################################
// Function Add_Uncut_Fine_Cell
//#####################################################################
template<class TV,int static_degree> void BASIS_INTEGRATION_UNIFORM_NEW<TV,static_degree>::
Add_Uncut_Fine_Cell(const TV_INT& cell,int block,int inside)
{
    for(int i=0;i<volume_blocks.m;i++)
        volume_blocks(i)->Add_Open_Subcell_Entries(cdi.Flatten(cell),block,inside);
}
//#####################################################################
// Function Cut_Elements
//#####################################################################
template<class TV,int static_degree> void BASIS_INTEGRATION_UNIFORM_NEW<TV,static_degree>::
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
// Function Add_Cut_Subcell
//#####################################################################
template<class TV,int static_degree> void BASIS_INTEGRATION_UNIFORM_NEW<TV,static_degree>::
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
            volume_inside+=Volume(projected_elements(i).x.X,dir);
        bool filled=volume_inside>subcell_range.Size()/2;
        Add_Uncut_Fine_Cell(cell,block,enclose_inside==filled);
        return;}

    STATIC_TENSOR<T,TV::m,static_degree+1> precomputed_integrals;
    RANGE<TV_INT> range(TV_INT(),TV_INT()+static_degree+1);
    for(RANGE_ITERATOR<TV::m> it(range);it.Valid();it.Next())
        if(volume_monomials_needed(it.index)){
            STATIC_POLYNOMIAL<T,TV::m,static_degree+1> monomial;
            monomial.Set_Term(it.index,1);
            monomial=monomial.Integrate(dir);
            T integral=0;
            for(int i=0;i<projected_elements.m;i++){
                const VECTOR<TV,TV::m>& V=projected_elements(i).x.X;
                integral+=monomial.Integrate_Over_Primitive(V)*T_FACE::Normal(V)(dir);}
            precomputed_integrals(it.index)+=integral;}

    Add_Uncut_Fine_Cell(cell,block,!enclose_inside);

    for(int i=0;i<volume_blocks.m;i++){
        VOLUME_BLOCK* vb=volume_blocks(i);
        for(int j=0;j<vb->overlap_polynomials.m;j++){
            typename VOLUME_BLOCK::OVERLAP_POLYNOMIAL& op=vb->overlap_polynomials(j);
            if(op.subcell&(1<<block)){
                T integral=Precomputed_Integral(precomputed_integrals,op.polynomial);
                int flat_index=cdi.Flatten(cell)+op.flat_index_offset;
                vb->Add_Entry(flat_index,op.flat_index_diff_ref,enclose_inside,integral);
                vb->Add_Entry(flat_index,op.flat_index_diff_ref,!enclose_inside,-integral);}}}

    STATIC_TENSOR<T,TV::m,static_degree+1> precomputed_interface_integrals[subcell_elements];
    bool has_element[subcell_elements]={};
    for(RANGE_ITERATOR<TV::m> it(range);it.Valid();it.Next())
        if(surface_monomials_needed(it.index)){
            STATIC_POLYNOMIAL<T,TV::m,static_degree> monomial;
            monomial.Set_Term(it.index,1);
            for(int i=0;i<interface_elements.m;i++){
                VECTOR<TV,TV::m> V=interface_elements(i).x.X;
                for(int j=0;j<TV::m;j++) V(j)=(V(j)*coarse_factor-TV(subcell_cell)-(T).5)*grid.dX;
                int e=interface_elements(i).y;
                has_element[e]=true;
                precomputed_interface_integrals[e](it.index)+=monomial.Integrate_Over_Primitive(V);}}
    
    for(int i=0;i<interface_blocks.m;i++){
        INTERFACE_BLOCK* ib=interface_blocks(i);
        int sign1=(ib->ignore_orientation||enclose_inside)?1:-1;
        int sign2=(ib->ignore_orientation)?1:-1;
        for(int j=0;j<ib->overlap_polynomials.m;j++){
            typename INTERFACE_BLOCK::OVERLAP_POLYNOMIAL& op=ib->overlap_polynomials(j);
            if(op.subcell&(1<<block))
                for(int k=0;k<subcell_elements;k++)
                    if(has_element[k]){
                        T integral=Precomputed_Integral(precomputed_interface_integrals[k],op.polynomial)*sign1;
                        ib->Add_Entry(element_base+k,op.flat_index_diff_ref,enclose_inside,integral);
                        ib->Add_Entry(element_base+k,op.flat_index_diff_ref,!enclose_inside,sign2*integral);}}}
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
