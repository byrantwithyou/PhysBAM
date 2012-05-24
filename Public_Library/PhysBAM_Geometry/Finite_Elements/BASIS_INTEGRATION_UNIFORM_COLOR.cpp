//#####################################################################
// Copyright 2012, Craig Schroeder, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <PhysBAM_Tools/Data_Structures/TRIPLE.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Matrices/MATRIX.h>
#include <PhysBAM_Geometry/Basic_Geometry/LINE_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry/PLANE.h>
#include <PhysBAM_Geometry/Basic_Geometry/SEGMENT_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <PhysBAM_Geometry/Finite_Elements/BASIS_INTEGRATION_UNIFORM_COLOR.h>
#include <PhysBAM_Geometry/Finite_Elements/CELL_DOMAIN_INTERFACE_COLOR.h>
#include <PhysBAM_Geometry/Grids_Uniform_Computations/MARCHING_CUBES_COLOR.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV,int static_degree> BASIS_INTEGRATION_UNIFORM_COLOR<TV,static_degree>::
BASIS_INTEGRATION_UNIFORM_COLOR(const GRID<TV>& grid_input,const GRID<TV>& phi_grid_input,const ARRAY<T,TV_INT>& phi_value_input,
    const ARRAY<int,TV_INT>& phi_color_input,CELL_DOMAIN_INTERFACE_COLOR<TV>& cdi_input)
    :grid(grid_input),phi_grid(phi_grid_input),phi_value(phi_value_input),phi_color(phi_color_input),cdi(cdi_input)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV,int static_degree> BASIS_INTEGRATION_UNIFORM_COLOR<TV,static_degree>::
~BASIS_INTEGRATION_UNIFORM_COLOR()
{
    volume_blocks.Delete_Pointers_And_Clean_Memory();
    surface_blocks.Delete_Pointers_And_Clean_Memory();
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
// Function Compute_Averaged_Orientation_Helper
//#####################################################################
template<class T,class T_FACE> static void
Compute_Averaged_Orientation_Helper(const VECTOR<ARRAY<TRIPLE<T_FACE,int,int> >,4>& surface,const HASHTABLE<VECTOR<int,2>,int>& ht_color_pairs,ARRAY<MATRIX<T,2> >& base_orientation)
{
}
//#####################################################################
// Function Compute_Averaged_Orientation_Helper
//#####################################################################
template<class T, class T_FACE> static void
Compute_Averaged_Orientation_Helper(const VECTOR<ARRAY<TRIPLE<T_FACE,int,int> >,8>& surface,const HASHTABLE<VECTOR<int,2>,int>& ht_color_pairs,ARRAY<MATRIX<T,3> >& base_orientation)
{
    ARRAY<VECTOR<T,3> > normal(base_orientation.m);
    VECTOR<T,3> tangent;
    for(int s=0;s<8;s++){
        const ARRAY<TRIPLE<T_FACE,int,int> >& subcell_surface=surface(s);
        for(int i=0;i<subcell_surface.m;i++){
            const TRIPLE<T_FACE,int,int>& surface_element=subcell_surface(i);
            if(surface_element.z<0) continue;
            int color_pair_index=-1;
            bool found=ht_color_pairs.Get(VECTOR<int,2>(surface_element.y,surface_element.z),color_pair_index);
            PHYSBAM_ASSERT(found);
            normal(color_pair_index)+=surface_element.x.Raw_Normal();}}

    for(int i=0;i<normal.m;i++){
        normal(i).Normalize();
        tangent=normal(i).Unit_Orthogonal_Vector();
        base_orientation(i)=MATRIX<T,3>(VECTOR<T,3>::Cross_Product(tangent,normal(i)),tangent,normal(i));}
}
//#####################################################################
// Function Compute_Entries
//#####################################################################
template<class TV,int static_degree> void BASIS_INTEGRATION_UNIFORM_COLOR<TV,static_degree>::
Compute_Entries(VECTOR<ARRAY<VECTOR_ND<T> >,TV::m>& f_surface)
{
    const VECTOR<TV_INT,(1<<TV::m)>& bits=GRID<TV>::Binary_Counts(TV_INT());
    MARCHING_CUBES_COLOR<TV>::Initialize_Case_Table();

    VECTOR<ARRAY<TRIPLE<T_FACE,int,int> >,(1<<TV::m)> surface;
    VECTOR<ARRAY<PAIR<T_FACE,int> >,(1<<TV::m)> sides;

    Compute_Open_Entries();

    for(UNIFORM_GRID_ITERATOR_CELL<TV> it(grid);it.Valid();it.Next()){
        TV_INT cell_base=it.index*2;
        int base_color=phi_color(cell_base);
        int cell_corners=0;
        bool material_cell=false;
        for(int b=0;b<(1<<TV::m);b++){
            int color=phi_color(cell_base+bits(b)*2);
            cell_corners|=(color!=base_color)<<b;
            material_cell|=(color>=0);}
        if(!material_cell) continue;
        if(cell_corners==0){Add_Uncut_Cell(it.index,base_color);continue;}

        VECTOR<bool,(1<<TV::m)> material_subcell;
        for(int s=0;s<(1<<TV::m);s++){
            sides(s).Remove_All();
            surface(s).Remove_All();
            TV_INT subcell_base=cell_base+bits(s);
            VECTOR<int,1<<TV::m> subcell_phi_colors;
            VECTOR<T,1<<TV::m> subcell_phi_values;
            for(int b=0;b<(1<<TV::m);b++){
                TV_INT subcell_vertex=subcell_base+bits(b);
                subcell_phi_colors(b)=phi_color(subcell_vertex);
                subcell_phi_values(b)=phi_value(subcell_vertex);
                material_subcell(s)|=subcell_phi_colors(b)>=0;}
            if(material_subcell(s)) MARCHING_CUBES_COLOR<TV>::Get_Elements_For_Cell(surface(s),sides(s),subcell_phi_colors,subcell_phi_values);}

        HASHTABLE<VECTOR<int,2>,int> ht_color_pairs;
        ARRAY<VECTOR<int,2> > color_pairs;
        
        for(int s=0;s<(1<<TV::m);s++){
            const ARRAY<TRIPLE<T_FACE,int,int> >& subcell_surface=surface(s);
            for(int i=0;i<subcell_surface.m;i++){
                const TRIPLE<T_FACE,int,int>& surface_element=subcell_surface(i);
                if(surface_element.z>=0){
                    VECTOR<int,2> color_pair(surface_element.y,surface_element.z);
                    if(!ht_color_pairs.Contains(color_pair))
                        ht_color_pairs.Insert(color_pair,color_pairs.Append(color_pair));}}}

        ARRAY<MATRIX<T,TV::m> > base_orientation(color_pairs.m);
        Compute_Averaged_Orientation_Helper(surface,ht_color_pairs,base_orientation);

        ARRAY<int> constraint_offsets(color_pairs.m);
        int full_constraints=0;
        int slip_constraints=0;
        for(int i=0;i<color_pairs.m;i++)
            if(color_pairs(i).x==-2||color_pairs(i).x>=0)
                constraint_offsets(i)=full_constraints++;
        for(int i=0;i<color_pairs.m;i++)
            if(color_pairs(i).x==-3)
                constraint_offsets(i)=slip_constraints+full_constraints++;

        cdi.Set_Flat_Base_And_Resize(full_constraints+slip_constraints,full_constraints,it.index);
        for(int i=0;i<surface_blocks.m;i++) surface_blocks(i)->Resize();

        for(int s=0;s<(1<<TV::m);s++)
            if(material_subcell(s)){
                if(!surface(s).m){
                    int color=phi_color(cell_base+bits(s));
                    assert(color>=0);
                    Add_Uncut_Fine_Cell(it.index,s,color);}
                else Add_Cut_Fine_Cell(it.index,s,TV(bits((1<<TV::m)-1-s)),surface(s),sides(s),
                    base_orientation,f_surface,constraint_offsets,ht_color_pairs);}
        cdi.Update_Constraint_Count();}
    cdi.Update_Total_Constraint_Count();
}
//#####################################################################
// Function Compute_Open_Entries
//#####################################################################
template<class TV,int static_degree> void BASIS_INTEGRATION_UNIFORM_COLOR<TV,static_degree>::
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
template<class TV,int static_degree> void BASIS_INTEGRATION_UNIFORM_COLOR<TV,static_degree>::
Add_Uncut_Cell(const TV_INT& cell,int color)
{
    for(int i=0;i<volume_blocks.m;i++)
        volume_blocks(i)->Add_Open_Entries(cdi.Flatten(cell),color);
}
//#####################################################################
// Function Add_Uncut_Fine_Cell
//#####################################################################
template<class TV,int static_degree> void BASIS_INTEGRATION_UNIFORM_COLOR<TV,static_degree>::
Add_Uncut_Fine_Cell(const TV_INT& cell,int subcell,int color)
{
    for(int i=0;i<volume_blocks.m;i++)
        volume_blocks(i)->Add_Open_Subcell_Entries(cdi.Flatten(cell),subcell,color);
}
//#####################################################################
// Function Compute_Consistent_Orientation_Helper
//#####################################################################
template<class T, class T_FACE> static void
Compute_Consistent_Orientation_Helper(const T_FACE& segment,MATRIX<T,2>& orientation,const MATRIX<T,2>& base_orientation)
{
    VECTOR<T,2> tangent=(segment.X(1)-segment.X(0)).Normalized();
    orientation=MATRIX<T,2>(tangent,tangent.Rotate_Clockwise_90());
}
//#####################################################################
// Function Compute_Consistent_Orientation_Helper
//#####################################################################
template<class T, class T_FACE> static void
Compute_Consistent_Orientation_Helper(const T_FACE& triangle,MATRIX<T,3>& orientation,const MATRIX<T,3>& base_orientation)
{
    typedef VECTOR<T,3> TV;
    typedef MATRIX<T,3> TM;

    TV b(0,0,1);
    TV n=base_orientation.Transpose_Times(triangle.Normal());
    TV u=TV::Cross_Product(b,n).Normalized();

    orientation=base_orientation*TM(n,u,TV::Cross_Product(n,u)).Times_Transpose(TM(b,u,TV::Cross_Product(b,u)));
}
//#####################################################################
// Function Add_Cut_Fine_Cell
//#####################################################################
template<class TV,int static_degree> void BASIS_INTEGRATION_UNIFORM_COLOR<TV,static_degree>::
Add_Cut_Fine_Cell(const TV_INT& cell,int subcell,const TV& subcell_offset,ARRAY<TRIPLE<T_FACE,int,int> >& surface,ARRAY<PAIR<T_FACE,int> >& sides,
    const ARRAY<MATRIX<T,TV::m> >& base_orientation,VECTOR<ARRAY<VECTOR_ND<T> >,TV::m>& f_surface,const ARRAY<int>& constraint_offsets,
    const HASHTABLE<VECTOR<int,2>,int>& ht_color_pairs)
{
    assert(sides.m);
    assert(surface.m);

    for(int i=0;i<surface.m;i++)
        for(int j=0;j<TV::m;j++)
            surface(i).x.X(j)=(surface(i).x.X(j)-subcell_offset)*((T).5*grid.dX);

    for(int i=0;i<sides.m;i++)
        for(int j=0;j<TV::m;j++)
            sides(i).x.X(j)=(sides(i).x.X(j)-subcell_offset)*((T).5*grid.dX);
            
    ARRAY<STATIC_TENSOR<T,TV::m,static_degree+1> > precomputed_volume_integrals(cdi.colors);
    RANGE<TV_INT> range(TV_INT(),TV_INT()+static_degree+1);
    ARRAY<T> integrals(cdi.colors);
    for(RANGE_ITERATOR<TV::m> it(range);it.Valid();it.Next())
        if(volume_monomials_needed(it.index)){
            STATIC_POLYNOMIAL<T,TV::m,static_degree+1> monomial;
            monomial.Set_Term(it.index,1);
            monomial=monomial.Integrate(TV::m-1);
            integrals.Fill(0);
            for(int i=0;i<sides.m;i++){
                const PAIR<T_FACE,int>& V=sides(i);
                if(V.y>=0) integrals(V.y)+=monomial.Quadrature_Over_Primitive(V.x.X)*T_FACE::Normal(V.x.X)(TV::m-1);}
            for(int i=0;i<surface.m;i++){
                const TRIPLE<T_FACE,int,int>& V=surface(i);
                T integral=monomial.Quadrature_Over_Primitive(V.x.X)*T_FACE::Normal(V.x.X)(TV::m-1);
                if(V.y>=0) integrals(V.y)-=integral;
                if(V.z>=0) integrals(V.z)+=integral;}
            for(int c=0;c<cdi.colors;c++) precomputed_volume_integrals(c)(it.index)+=integrals(c);}

    for(int i=0;i<volume_blocks.m;i++){
        VOLUME_BLOCK* vb=volume_blocks(i);
        for(int j=0;j<vb->overlap_polynomials.m;j++){
            typename VOLUME_BLOCK::OVERLAP_POLYNOMIAL& op=vb->overlap_polynomials(j);
            if(op.subcell&(1<<subcell))
                for(int c=0;c<cdi.colors;c++){
                    T integral=Precomputed_Integral(precomputed_volume_integrals(c),op.polynomial);
                    int flat_index=cdi.Flatten(cell)+op.flat_index_offset;
                    vb->Add_Entry(flat_index,op.flat_index_diff_ref,c,integral);}}}

    ARRAY<MATRIX<T,TV::m> > orientations;
    orientations.Resize(surface.m);
    for(int i=0;i<surface.m;i++){
        const TRIPLE<T_FACE,int,int>& surface_element=surface(i);
        if(surface_element.z<0) continue;
        int color_pair_index=-1;
        bool found=ht_color_pairs.Get(VECTOR<int,2>(surface_element.y,surface_element.z),color_pair_index);
        PHYSBAM_ASSERT(found);
        Compute_Consistent_Orientation_Helper(surface(i).x,orientations(i),base_orientation(color_pair_index));}

    ARRAY<STATIC_TENSOR<T,TV::m,static_degree+1> > precomputed_surface_integrals(surface.m);
    for(RANGE_ITERATOR<TV::m> it(range);it.Valid();it.Next())
        if(surface_monomials_needed(it.index)){
            STATIC_POLYNOMIAL<T,TV::m,static_degree> monomial;
            monomial.Set_Term(it.index,1);
            for(int k=0;k<surface.m;k++) precomputed_surface_integrals(k)(it.index)=monomial.Quadrature_Over_Primitive(surface(k).x.X);}

    for(int i=0;i<surface_blocks.m;i++){
        SURFACE_BLOCK* sb=surface_blocks(i);
        for(int j=0;j<sb->overlap_polynomials.m;j++){
            typename SURFACE_BLOCK::OVERLAP_POLYNOMIAL& op=sb->overlap_polynomials(j);
            if(op.subcell&(1<<subcell))
                for(int k=0;k<surface.m;k++){
                    const TRIPLE<T_FACE,int,int>& surface_element=surface(k);
                    if(surface_element.z<0) continue;
                    int color_pair_index=-1;
                    bool found=ht_color_pairs.Get(VECTOR<int,2>(surface_element.y,surface_element.z),color_pair_index);
                    PHYSBAM_ASSERT(found);
                    T integral=Precomputed_Integral(precomputed_surface_integrals(k),op.polynomial);
                    int constraint_offset=constraint_offsets(color_pair_index);
                        
                        if(surface_element.y==-2||surface_element.y>=0)
                            for(int orientation=0;orientation<TV::m-1;orientation++){
                                T value=integral*orientations(k)(sb->axis,orientation);
                                if(surface_element.y>=0) sb->Add_Entry(cdi.constraint_base_tangent+constraint_offset,orientation,op.flat_index_diff_ref,surface_element.y,value);
                                if(surface_element.z>=0) sb->Add_Entry(cdi.constraint_base_tangent+constraint_offset,orientation,op.flat_index_diff_ref,surface_element.z,-value);}

                        if(surface_element.y!=-1){
                            T value=integral*orientations(k)(sb->axis,TV::m-1);
                            if(surface_element.y>=0) sb->Add_Entry(cdi.constraint_base_normal+constraint_offset,TV::m-1,op.flat_index_diff_ref,surface_element.y,value);
                            if(surface_element.z>=0) sb->Add_Entry(cdi.constraint_base_normal+constraint_offset,TV::m-1,op.flat_index_diff_ref,surface_element.z,-value);}

                        if(surface_element.y>=0||surface_element.y==-1){
                            int flat_index=cdi.Flatten(cell)+sb->Flat_Diff(op.flat_index_diff_ref);
                            T value=integral*sb->abc->f_surface(surface_element.x.Center()+grid.Center(cell),surface_element.y,surface_element.z)(sb->axis);
                            if(surface_element.y>=0) value*=-0.5;
                            if(surface_element.y>=0) f_surface(sb->axis)(surface_element.y)(flat_index)+=value;
                            if(surface_element.z>=0) f_surface(sb->axis)(surface_element.z)(flat_index)+=value;}}}}
}
//#####################################################################
// Function Add_Volume_Block
//#####################################################################
template<class TV,int static_degree> template<int d0,int d1> void BASIS_INTEGRATION_UNIFORM_COLOR<TV,static_degree>::
Add_Volume_Block(SYSTEM_VOLUME_BLOCK_HELPER_COLOR<TV>& helper,const BASIS_STENCIL_UNIFORM<TV,d0>& s0,
    const BASIS_STENCIL_UNIFORM<TV,d1>& s1,const ARRAY<T>& scale)
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
// Function Add_Surface_Block
//#####################################################################
template<class TV,int static_degree> template<int d> void BASIS_INTEGRATION_UNIFORM_COLOR<TV,static_degree>::
Add_Surface_Block(SYSTEM_SURFACE_BLOCK_HELPER_COLOR<TV>& helper,const BASIS_STENCIL_UNIFORM<TV,d>& s,
    ANALYTIC_BOUNDARY_CONDITIONS_COLOR<TV>* abc,int axis,T scale)
{
    SURFACE_BLOCK* sb=new SURFACE_BLOCK;
    sb->Initialize(helper,s,abc,axis,scale);
    surface_blocks.Append(sb);
        
    for(int i=0;i<sb->overlap_polynomials.m;i++){
        RANGE<TV_INT> range(TV_INT(),sb->overlap_polynomials(i).polynomial.size+1);
        for(RANGE_ITERATOR<TV::m> it(range);it.Valid();it.Next())
            if(sb->overlap_polynomials(i).polynomial.terms(it.index))
                surface_monomials_needed(it.index)=true;}
}
template class BASIS_INTEGRATION_UNIFORM_COLOR<VECTOR<float,3>,2>;
template class BASIS_INTEGRATION_UNIFORM_COLOR<VECTOR<float,2>,2>;
template void BASIS_INTEGRATION_UNIFORM_COLOR<VECTOR<float,3>,2>::Add_Volume_Block<0,1>(SYSTEM_VOLUME_BLOCK_HELPER_COLOR<VECTOR<float,3> >&,
    BASIS_STENCIL_UNIFORM<VECTOR<float,3>,0> const&,BASIS_STENCIL_UNIFORM<VECTOR<float,3>,1> const&,ARRAY<float> const&);
template void BASIS_INTEGRATION_UNIFORM_COLOR<VECTOR<float,3>,2>::Add_Volume_Block<1,1>(SYSTEM_VOLUME_BLOCK_HELPER_COLOR<VECTOR<float,3> >&,
    BASIS_STENCIL_UNIFORM<VECTOR<float,3>,1> const&,BASIS_STENCIL_UNIFORM<VECTOR<float,3>,1> const&,ARRAY<float> const&);
template void BASIS_INTEGRATION_UNIFORM_COLOR<VECTOR<float,3>,2>::Add_Surface_Block<1>(SYSTEM_SURFACE_BLOCK_HELPER_COLOR<VECTOR<float,3> >&,
    BASIS_STENCIL_UNIFORM<VECTOR<float,3>,1> const&,ANALYTIC_BOUNDARY_CONDITIONS_COLOR<VECTOR<float,3> >*,int,float);
template void BASIS_INTEGRATION_UNIFORM_COLOR<VECTOR<float,2>,2>::Add_Volume_Block<0,1>(SYSTEM_VOLUME_BLOCK_HELPER_COLOR<VECTOR<float,2> >&,
    BASIS_STENCIL_UNIFORM<VECTOR<float,2>,0> const&,BASIS_STENCIL_UNIFORM<VECTOR<float,2>,1> const&,ARRAY<float> const&);
template void BASIS_INTEGRATION_UNIFORM_COLOR<VECTOR<float,2>,2>::Add_Volume_Block<1,1>(SYSTEM_VOLUME_BLOCK_HELPER_COLOR<VECTOR<float,2> >&,
    BASIS_STENCIL_UNIFORM<VECTOR<float,2>,1> const&,BASIS_STENCIL_UNIFORM<VECTOR<float,2>,1> const&,ARRAY<float> const&);
template void BASIS_INTEGRATION_UNIFORM_COLOR<VECTOR<float,2>,2>::Add_Surface_Block<1>(SYSTEM_SURFACE_BLOCK_HELPER_COLOR<VECTOR<float,2> >&,
    BASIS_STENCIL_UNIFORM<VECTOR<float,2>,1> const&,ANALYTIC_BOUNDARY_CONDITIONS_COLOR<VECTOR<float,2> >*,int,float);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class BASIS_INTEGRATION_UNIFORM_COLOR<VECTOR<double,3>,2>;
template class BASIS_INTEGRATION_UNIFORM_COLOR<VECTOR<double,2>,2>;
template void BASIS_INTEGRATION_UNIFORM_COLOR<VECTOR<double,3>,2>::Add_Volume_Block<0,1>(SYSTEM_VOLUME_BLOCK_HELPER_COLOR<VECTOR<double,3> >&,
    BASIS_STENCIL_UNIFORM<VECTOR<double,3>,0> const&,BASIS_STENCIL_UNIFORM<VECTOR<double,3>,1> const&,ARRAY<double> const&);
template void BASIS_INTEGRATION_UNIFORM_COLOR<VECTOR<double,3>,2>::Add_Volume_Block<1,1>(SYSTEM_VOLUME_BLOCK_HELPER_COLOR<VECTOR<double,3> >&,
    BASIS_STENCIL_UNIFORM<VECTOR<double,3>,1> const&,BASIS_STENCIL_UNIFORM<VECTOR<double,3>,1> const&,ARRAY<double> const&);
template void BASIS_INTEGRATION_UNIFORM_COLOR<VECTOR<double,3>,2>::Add_Surface_Block<1>(SYSTEM_SURFACE_BLOCK_HELPER_COLOR<VECTOR<double,3> >&,
    BASIS_STENCIL_UNIFORM<VECTOR<double,3>,1> const&,ANALYTIC_BOUNDARY_CONDITIONS_COLOR<VECTOR<double,3> >*,int,double);
template void BASIS_INTEGRATION_UNIFORM_COLOR<VECTOR<double,2>,2>::Add_Volume_Block<0,1>(SYSTEM_VOLUME_BLOCK_HELPER_COLOR<VECTOR<double,2> >&,
    BASIS_STENCIL_UNIFORM<VECTOR<double,2>,0> const&,BASIS_STENCIL_UNIFORM<VECTOR<double,2>,1> const&,ARRAY<double> const&);
template void BASIS_INTEGRATION_UNIFORM_COLOR<VECTOR<double,2>,2>::Add_Volume_Block<1,1>(SYSTEM_VOLUME_BLOCK_HELPER_COLOR<VECTOR<double,2> >&,
    BASIS_STENCIL_UNIFORM<VECTOR<double,2>,1> const&,BASIS_STENCIL_UNIFORM<VECTOR<double,2>,1> const&,ARRAY<double> const&);
template void BASIS_INTEGRATION_UNIFORM_COLOR<VECTOR<double,2>,2>::Add_Surface_Block<1>(SYSTEM_SURFACE_BLOCK_HELPER_COLOR<VECTOR<double,2> >&,
    BASIS_STENCIL_UNIFORM<VECTOR<double,2>,1> const&,ANALYTIC_BOUNDARY_CONDITIONS_COLOR<VECTOR<double,2> >*,int,double);
#endif
