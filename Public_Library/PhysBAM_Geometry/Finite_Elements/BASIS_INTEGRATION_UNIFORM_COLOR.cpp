//#####################################################################
// Copyright 2012, Craig Schroeder, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Matrices/MATRIX.h>
#include <PhysBAM_Geometry/Finite_Elements/BASIS_INTEGRATION_UNIFORM_COLOR.h>
#include <PhysBAM_Geometry/Finite_Elements/CELL_DOMAIN_INTERFACE_COLOR.h>
#include <PhysBAM_Geometry/Grids_Uniform_Computations/MARCHING_CUBES_COLOR.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>

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
    surface_blocks_scalar.Delete_Pointers_And_Clean_Memory();
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
template<class T,class CELL_ELEMENTS> static void
Compute_Averaged_Orientation_Helper(const VECTOR<CELL_ELEMENTS,4>& cell_elements,const HASHTABLE<VECTOR<int,2>,int>& ht_color_pairs,ARRAY<MATRIX<T,2> >& base_orientation)
{
}
//#####################################################################
// Function Compute_Averaged_Orientation_Helper
//#####################################################################
template<class T, class CELL_ELEMENTS> static void
Compute_Averaged_Orientation_Helper(const VECTOR<CELL_ELEMENTS,8>& cell_elements,const HASHTABLE<VECTOR<int,2>,int>& ht_color_pairs,ARRAY<MATRIX<T,3> >& base_orientation)
{
    typedef VECTOR<T,3> TV;
    typedef typename CELL_DOMAIN_INTERFACE_COLOR<TV>::INTERFACE_ELEMENT INTERFACE_ELEMENT;

    ARRAY<TV > normal(base_orientation.m);
    TV tangent;
    for(int s=0;s<8;s++){
        const ARRAY<INTERFACE_ELEMENT>& subcell_interface=cell_elements(s).interface;
        for(int i=0;i<subcell_interface.m;i++){
            const INTERFACE_ELEMENT& V=subcell_interface(i);
            if(V.color_pair.y<0) continue;
            int color_pair_index=-1;
            bool found=ht_color_pairs.Get(VECTOR<int,2>(V.color_pair.x,V.color_pair.y),color_pair_index);
            PHYSBAM_ASSERT(found);
            normal(color_pair_index)+=V.face.Raw_Normal();}}

    for(int i=0;i<normal.m;i++){
        normal(i).Normalize();
        tangent=normal(i).Unit_Orthogonal_Vector();
        base_orientation(i)=MATRIX<T,3>(TV::Cross_Product(tangent,normal(i)),tangent,normal(i));}
}
//#####################################################################
// Function Compute_Entries
//#####################################################################
template<class TV,int static_degree> void BASIS_INTEGRATION_UNIFORM_COLOR<TV,static_degree>::
Compute_Entries()
{
    const VECTOR<TV_INT,(1<<TV::m)>& bits=GRID<TV>::Binary_Counts(TV_INT());

    Compute_Open_Entries();

    for(UNIFORM_GRID_ITERATOR_CELL<TV> it(grid);it.Valid();it.Next()){
        VECTOR<CELL_ELEMENTS,(1<<TV::m)> cell_elements;
        TV_INT cell_base=it.index*2;
        TV cell_center(grid.Center(it.index));
        for(int s=0;s<(1<<TV::m);s++)
            if(const CELL_ELEMENTS* pce=cdi.index_to_cell_elements.Get_Pointer(cell_base+bits(s))){
                CELL_ELEMENTS& ce=(cell_elements(s)=*pce);
                for(int i=0;i<ce.interface.m;i++)
                    ce.interface(i).face.X-=cell_center;
                for(int i=0;i<ce.boundary.m;i++)
                    ce.boundary(i).face.X-=cell_center;}

        HASHTABLE<VECTOR<int,2>,int> ht_color_pairs;
        ARRAY<VECTOR<int,2> > color_pairs;
        
        for(int s=0;s<(1<<TV::m);s++){
            const ARRAY<INTERFACE_ELEMENT>& subcell_interface=cell_elements(s).interface;
            for(int i=0;i<subcell_interface.m;i++){
                const INTERFACE_ELEMENT& V=subcell_interface(i);
                if(V.color_pair.y>=0){
                    VECTOR<int,2> color_pair(V.color_pair.x,V.color_pair.y);
                    if(!ht_color_pairs.Contains(color_pair))
                        ht_color_pairs.Insert(color_pair,color_pairs.Append(color_pair));}}}

        ARRAY<int> constraint_offsets(color_pairs.m);
        int full_constraints=0;
        int slip_constraints=0;
        for(int i=0;i<color_pairs.m;i++)
            if(color_pairs(i).x==BC::DIRICHLET||color_pairs(i).x>=0)
                constraint_offsets(i)=full_constraints++;
        for(int i=0;i<color_pairs.m;i++)
            if(color_pairs(i).x==BC::SLIP)
                constraint_offsets(i)=slip_constraints+full_constraints++;
        
        all_constraint_color_pairs.Append(color_pairs);
        all_constraint_offsets.Append(constraint_offsets);

        ARRAY<MATRIX<T,TV::m> > base_orientation;
        if(surface_blocks.m){
            base_orientation.Resize(color_pairs.m);
            Compute_Averaged_Orientation_Helper(cell_elements,ht_color_pairs,base_orientation);
            cdi.Set_Flat_Base_And_Resize(full_constraints+slip_constraints,full_constraints,it.index);
            for(int i=0;i<surface_blocks.m;i++) surface_blocks(i)->Resize();}
        if(surface_blocks_scalar.m){
            cdi.Set_Flat_Base_And_Resize_Scalar(full_constraints,it.index);
            for(int i=0;i<surface_blocks_scalar.m;i++) surface_blocks_scalar(i)->Resize();
            assert(slip_constraints==0 && "Do you really want slip constraints for a scalar variable?");}

        for(int s=0;s<(1<<TV::m);s++){
            if(!cell_elements(s).interface.m){
                int color=phi_color(cell_base+bits(s));
                if(color>=0)
                    Add_Uncut_Fine_Cell(it.index,s,color);}
            else Add_Cut_Fine_Cell(it.index,s,TV(bits((1<<TV::m)-1-s)),cell_elements(s),
                base_orientation,constraint_offsets,ht_color_pairs);}
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
Add_Cut_Fine_Cell(const TV_INT& cell,int subcell,const TV& subcell_offset,const CELL_ELEMENTS& cell_elements,
    const ARRAY<MATRIX<T,TV::m> >& base_orientation,const ARRAY<int>& constraint_offsets,const HASHTABLE<VECTOR<int,2>,int>& ht_color_pairs)
{
    assert(cell_elements.boundary.m);
    assert(cell_elements.interface.m);

    TV cell_center(grid.Center(cell));
    
    ARRAY<STATIC_TENSOR<T,TV::m,static_degree+1> > precomputed_volume_integrals(cdi.colors);
    RANGE<TV_INT> range(TV_INT(),TV_INT()+static_degree+1);
    ARRAY<T> integrals(cdi.colors);
    for(RANGE_ITERATOR<TV::m> it(range);it.Valid();it.Next())
        if(volume_monomials_needed(it.index)){
            STATIC_POLYNOMIAL<T,TV::m,static_degree+1> monomial;
            monomial.Set_Term(it.index,1);
            monomial=monomial.Integrate(TV::m-1);
            integrals.Fill(0);
            for(int i=0;i<cell_elements.boundary.m;i++){
                const BOUNDARY_ELEMENT& V=cell_elements.boundary(i);
                if(V.color>=0) integrals(V.color)+=monomial.Quadrature_Over_Primitive(V.face.X)*T_FACE::Normal(V.face.X)(TV::m-1);}
            for(int i=0;i<cell_elements.interface.m;i++){
                const INTERFACE_ELEMENT& V=cell_elements.interface(i);
                if(V.color_pair.y<0) continue;
                T integral=monomial.Quadrature_Over_Primitive(V.face.X)*T_FACE::Normal(V.face.X)(TV::m-1);
                if(V.color_pair.x>=0) integrals(V.color_pair.x)-=integral;
                if("Alexey was here") integrals(V.color_pair.y)+=integral;}
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

    ARRAY<STATIC_TENSOR<T,TV::m,static_degree+1> > precomputed_surface_integrals(cell_elements.interface.m);
    for(RANGE_ITERATOR<TV::m> it(range);it.Valid();it.Next())
        if(surface_monomials_needed(it.index)){
            STATIC_POLYNOMIAL<T,TV::m,static_degree> monomial;
            monomial.Set_Term(it.index,1);
            for(int k=0;k<cell_elements.interface.m;k++) precomputed_surface_integrals(k)(it.index)=monomial.Quadrature_Over_Primitive(cell_elements.interface(k).face.X);}

    if(surface_blocks.m){
        ARRAY<MATRIX<T,TV::m> > orientations(cell_elements.interface.m);
        for(int i=0;i<cell_elements.interface.m;i++){
            const INTERFACE_ELEMENT& V=cell_elements.interface(i);
            if(V.color_pair.y<0) continue;
            int color_pair_index=-1;
            bool found=ht_color_pairs.Get(VECTOR<int,2>(V.color_pair.x,V.color_pair.y),color_pair_index);
            PHYSBAM_ASSERT(found);
            Compute_Consistent_Orientation_Helper(cell_elements.interface(i).face,orientations(i),base_orientation(color_pair_index));}

        for(int i=0;i<surface_blocks.m;i++){
            SURFACE_BLOCK* sb=surface_blocks(i);
            for(int j=0;j<sb->overlap_polynomials.m;j++){
                typename SURFACE_BLOCK::OVERLAP_POLYNOMIAL& op=sb->overlap_polynomials(j);
                if(op.subcell&(1<<subcell))
                    for(int k=0;k<cell_elements.interface.m;k++){
                        const INTERFACE_ELEMENT& V=cell_elements.interface(k);
                        if(V.color_pair.y<0) continue;
                        int color_pair_index=-1;
                        bool found=ht_color_pairs.Get(VECTOR<int,2>(V.color_pair.x,V.color_pair.y),color_pair_index);
                        PHYSBAM_ASSERT(found);
                        const int constraint_offset=constraint_offsets(color_pair_index);
                        const T integral=Precomputed_Integral(precomputed_surface_integrals(k),op.polynomial);
                        int flat_index=cdi.Flatten(cell)+sb->Flat_Diff(op.flat_index_diff_ref);

                        cdi.nc_present|=(V.color_pair.x==BC::NEUMANN);
                        cdi.dc_present|=(V.color_pair.x==BC::DIRICHLET);
                        cdi.sc_present|=(V.color_pair.x==BC::SLIP);

                        if(V.color_pair.x!=BC::SLIP && V.color_pair.x!=BC::NEUMANN)
                            for(int orientation=0;orientation<TV::m-1;orientation++){
                                T value=integral*orientations(k)(sb->axis,orientation);
                                if(V.color_pair.x>=0) sb->Add_Entry(cdi.constraint_base_t+constraint_offset,orientation,op.flat_index_diff_ref,V.color_pair.x,value);
                                if("Alexey was here") sb->Add_Entry(cdi.constraint_base_t+constraint_offset,orientation,op.flat_index_diff_ref,V.color_pair.y,-value);}
                        
                        if(V.color_pair.x!=BC::NEUMANN){
                            T value=integral*orientations(k)(sb->axis,TV::m-1);
                            if(V.color_pair.x>=0) sb->Add_Entry(cdi.constraint_base_n+constraint_offset,TV::m-1,op.flat_index_diff_ref,V.color_pair.x,value);
                            if("Alexey was here") sb->Add_Entry(cdi.constraint_base_n+constraint_offset,TV::m-1,op.flat_index_diff_ref,V.color_pair.y,-value);}
                        
                        TV X=V.face.Center()+grid.Center(cell);
                        if(V.color_pair.x>=0){
                            T value=-integral*sb->bc->j_surface(X,V.color_pair.x,V.color_pair.y)(sb->axis);
                            value*=(T)0.5;
                            (*sb->rhs)(V.color_pair.x)(flat_index)+=value;
                            (*sb->rhs)(V.color_pair.y)(flat_index)+=value;

                            if(sb->bc->use_discontinuous_velocity)
                                if(sb->axis==0){ // This code should not be repeated for each block
                                    TV value=-integral*orientations(k).Transpose_Times(sb->bc->u_jump(X,V.color_pair.x,V.color_pair.y));
                                    for(int d=0;d<TV::m;d++)
                                        sb->Add_Constraint_Rhs_Entry(*cdi.constraint_base(d)+constraint_offset,d,V.color_pair.y,value(d));}}
                        else if(V.color_pair.x==BC::NEUMANN){
                            T value=integral*sb->bc->n_surface(X,V.color_pair.x,V.color_pair.y)(sb->axis);
                            (*sb->rhs)(V.color_pair.y)(flat_index)+=value;}
                        else if(V.color_pair.x==BC::DIRICHLET){
                            if(sb->axis==0){ // This code should not be repeated for each block
                                TV value=-integral*orientations(k).Transpose_Times(sb->bc->d_surface(X,V.color_pair.x,V.color_pair.y));
                                for(int d=0;d<TV::m;d++)
                                    sb->Add_Constraint_Rhs_Entry(*cdi.constraint_base(d)+constraint_offset,d,V.color_pair.y,value(d));}}
                        else{
                            TV N=orientations(k).Column(TV::m-1);
                            T n_value=integral*sb->bc->n_surface(X,V.color_pair.x,V.color_pair.y).Projected_Orthogonal_To_Unit_Direction(N)(sb->axis);
                            (*sb->rhs)(V.color_pair.y)(flat_index)+=n_value;
                            if(sb->axis==0){ // This code should not be repeated for each block
                                T d_value=-integral*sb->bc->d_surface(X,V.color_pair.x,V.color_pair.y).Dot(N);
                                sb->Add_Constraint_Rhs_Entry(cdi.constraint_base_n+constraint_offset,TV::m-1,V.color_pair.y,d_value);}}}}}}

    if(surface_blocks_scalar.m){
        for(int i=0;i<surface_blocks_scalar.m;i++){
            SURFACE_BLOCK_SCALAR* sbs=surface_blocks_scalar(i);
            for(int j=0;j<sbs->overlap_polynomials.m;j++){
                typename SURFACE_BLOCK_SCALAR::OVERLAP_POLYNOMIAL& op=sbs->overlap_polynomials(j);
                if(op.subcell&(1<<subcell))
                    for(int k=0;k<cell_elements.interface.m;k++){
                        const INTERFACE_ELEMENT& V=cell_elements.interface(k);
                        if(V.color_pair.y<0) continue;
                        int color_pair_index=-1;
                        bool found=ht_color_pairs.Get(VECTOR<int,2>(V.color_pair.x,V.color_pair.y),color_pair_index);
                        PHYSBAM_ASSERT(found);
                        const int constraint_offset=constraint_offsets(color_pair_index);
                        const T integral=Precomputed_Integral(precomputed_surface_integrals(k),op.polynomial);
                        PHYSBAM_ASSERT(V.color_pair.y!=BC::SLIP && "Slip constraint is not allowed for scalar varialbes!");

                        cdi.nc_present|=(V.color_pair.x==BC::NEUMANN);
                        cdi.dc_present|=(V.color_pair.x==BC::DIRICHLET);
                        
                        if(V.color_pair.x!=BC::NEUMANN){
                            if(V.color_pair.x>=0) sbs->Add_Entry(cdi.constraint_base_scalar+constraint_offset,op.flat_index_diff_ref,V.color_pair.x,integral);
                            if("Alexey was here") sbs->Add_Entry(cdi.constraint_base_scalar+constraint_offset,op.flat_index_diff_ref,V.color_pair.y,-integral);}

                        TV X=V.face.Center()+cell_center;
                        int flat_index=cdi.Flatten(cell)+sbs->Flat_Diff(op.flat_index_diff_ref);
                        if(V.color_pair.x>=0){
                            T value=-integral*sbs->bc->j_surface(X,V.color_pair.x,V.color_pair.y);
                            value*=(T)0.5;
                            (*sbs->rhs)(V.color_pair.x)(flat_index)+=value;
                            (*sbs->rhs)(V.color_pair.y)(flat_index)+=value;
                            
                            if(sbs->bc->use_discontinuous_scalar_field){
                                T value=-integral*sbs->bc->u_jump(X,V.color_pair.x,V.color_pair.y);
                                sbs->Add_Constraint_Rhs_Entry(cdi.constraint_base_scalar+constraint_offset,V.color_pair.y,value);}}
                        else if(V.color_pair.x==BC::NEUMANN){
                            T value=-integral*sbs->bc->n_surface(X,V.color_pair.x,V.color_pair.y);
                            (*sbs->rhs)(V.color_pair.y)(flat_index)+=value;}
                        else{
                            T value=-integral*sbs->bc->d_surface(X,V.color_pair.x,V.color_pair.y);
                            sbs->Add_Constraint_Rhs_Entry(cdi.constraint_base_scalar+constraint_offset,V.color_pair.y,value);}}}}}
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
    BOUNDARY_CONDITIONS_COLOR<TV>* bc,ARRAY<ARRAY<T> >& f_surface,int axis,T scale)
{
    SURFACE_BLOCK* sb=new SURFACE_BLOCK;
    sb->Initialize(helper,s,bc,f_surface,axis,scale);
    surface_blocks.Append(sb);
        
    for(int i=0;i<sb->overlap_polynomials.m;i++){
        RANGE<TV_INT> range(TV_INT(),sb->overlap_polynomials(i).polynomial.size+1);
        for(RANGE_ITERATOR<TV::m> it(range);it.Valid();it.Next())
            if(sb->overlap_polynomials(i).polynomial.terms(it.index))
                surface_monomials_needed(it.index)=true;}
}
//#####################################################################
// Function Add_Surface_Block_Scalar
//#####################################################################
template<class TV,int static_degree> template<int d> void BASIS_INTEGRATION_UNIFORM_COLOR<TV,static_degree>::
Add_Surface_Block_Scalar(SYSTEM_SURFACE_BLOCK_SCALAR_HELPER_COLOR<TV>& helper,const BASIS_STENCIL_UNIFORM<TV,d>& s,
    BOUNDARY_CONDITIONS_SCALAR_COLOR<TV>* bc,ARRAY<ARRAY<T> >& f_surface,T scale)
{
    SURFACE_BLOCK_SCALAR* sbs=new SURFACE_BLOCK_SCALAR;
    sbs->Initialize(helper,s,bc,f_surface,scale);
    surface_blocks_scalar.Append(sbs);
        
    for(int i=0;i<sbs->overlap_polynomials.m;i++){
        RANGE<TV_INT> range(TV_INT(),sbs->overlap_polynomials(i).polynomial.size+1);
        for(RANGE_ITERATOR<TV::m> it(range);it.Valid();it.Next())
            if(sbs->overlap_polynomials(i).polynomial.terms(it.index))
                surface_monomials_needed(it.index)=true;}
}
namespace PhysBAM{
template class BASIS_INTEGRATION_UNIFORM_COLOR<VECTOR<float,3>,2>;
template class BASIS_INTEGRATION_UNIFORM_COLOR<VECTOR<float,2>,2>;
template void BASIS_INTEGRATION_UNIFORM_COLOR<VECTOR<float,3>,2>::Add_Volume_Block<0,1>(SYSTEM_VOLUME_BLOCK_HELPER_COLOR<VECTOR<float,3> >&,
    BASIS_STENCIL_UNIFORM<VECTOR<float,3>,0> const&,BASIS_STENCIL_UNIFORM<VECTOR<float,3>,1> const&,ARRAY<float> const&);
template void BASIS_INTEGRATION_UNIFORM_COLOR<VECTOR<float,3>,2>::Add_Volume_Block<1,1>(SYSTEM_VOLUME_BLOCK_HELPER_COLOR<VECTOR<float,3> >&,
    BASIS_STENCIL_UNIFORM<VECTOR<float,3>,1> const&,BASIS_STENCIL_UNIFORM<VECTOR<float,3>,1> const&,ARRAY<float> const&);
template void BASIS_INTEGRATION_UNIFORM_COLOR<VECTOR<float,3>,2>::Add_Surface_Block<1>(SYSTEM_SURFACE_BLOCK_HELPER_COLOR<VECTOR<float,3> >&,
    BASIS_STENCIL_UNIFORM<VECTOR<float,3>,1> const&,BOUNDARY_CONDITIONS_COLOR<VECTOR<float,3> >*,ARRAY<ARRAY<float> >&,int,float);
template void BASIS_INTEGRATION_UNIFORM_COLOR<VECTOR<float,3>,2>::Add_Surface_Block_Scalar<1>(SYSTEM_SURFACE_BLOCK_SCALAR_HELPER_COLOR<VECTOR<float,3> >&,
    BASIS_STENCIL_UNIFORM<VECTOR<float,3>,1> const&,BOUNDARY_CONDITIONS_SCALAR_COLOR<VECTOR<float,3> >*,ARRAY<ARRAY<float> >&,float);
template void BASIS_INTEGRATION_UNIFORM_COLOR<VECTOR<float,2>,2>::Add_Volume_Block<0,1>(SYSTEM_VOLUME_BLOCK_HELPER_COLOR<VECTOR<float,2> >&,
    BASIS_STENCIL_UNIFORM<VECTOR<float,2>,0> const&,BASIS_STENCIL_UNIFORM<VECTOR<float,2>,1> const&,ARRAY<float> const&);
template void BASIS_INTEGRATION_UNIFORM_COLOR<VECTOR<float,2>,2>::Add_Volume_Block<1,1>(SYSTEM_VOLUME_BLOCK_HELPER_COLOR<VECTOR<float,2> >&,
    BASIS_STENCIL_UNIFORM<VECTOR<float,2>,1> const&,BASIS_STENCIL_UNIFORM<VECTOR<float,2>,1> const&,ARRAY<float> const&);
template void BASIS_INTEGRATION_UNIFORM_COLOR<VECTOR<float,2>,2>::Add_Surface_Block<1>(SYSTEM_SURFACE_BLOCK_HELPER_COLOR<VECTOR<float,2> >&,
    BASIS_STENCIL_UNIFORM<VECTOR<float,2>,1> const&,BOUNDARY_CONDITIONS_COLOR<VECTOR<float,2> >*,ARRAY<ARRAY<float> >&,int,float);
template void BASIS_INTEGRATION_UNIFORM_COLOR<VECTOR<float,2>,2>::Add_Surface_Block_Scalar<1>(SYSTEM_SURFACE_BLOCK_SCALAR_HELPER_COLOR<VECTOR<float,2> >&,
    BASIS_STENCIL_UNIFORM<VECTOR<float,2>,1> const&,BOUNDARY_CONDITIONS_SCALAR_COLOR<VECTOR<float,2> >*,ARRAY<ARRAY<float> >&,float);
template class BASIS_INTEGRATION_UNIFORM_COLOR<VECTOR<double,3>,2>;
template class BASIS_INTEGRATION_UNIFORM_COLOR<VECTOR<double,2>,2>;
template void BASIS_INTEGRATION_UNIFORM_COLOR<VECTOR<double,3>,2>::Add_Volume_Block<0,1>(SYSTEM_VOLUME_BLOCK_HELPER_COLOR<VECTOR<double,3> >&,
    BASIS_STENCIL_UNIFORM<VECTOR<double,3>,0> const&,BASIS_STENCIL_UNIFORM<VECTOR<double,3>,1> const&,ARRAY<double> const&);
template void BASIS_INTEGRATION_UNIFORM_COLOR<VECTOR<double,3>,2>::Add_Volume_Block<1,1>(SYSTEM_VOLUME_BLOCK_HELPER_COLOR<VECTOR<double,3> >&,
    BASIS_STENCIL_UNIFORM<VECTOR<double,3>,1> const&,BASIS_STENCIL_UNIFORM<VECTOR<double,3>,1> const&,ARRAY<double> const&);
template void BASIS_INTEGRATION_UNIFORM_COLOR<VECTOR<double,3>,2>::Add_Surface_Block<1>(SYSTEM_SURFACE_BLOCK_HELPER_COLOR<VECTOR<double,3> >&,
    BASIS_STENCIL_UNIFORM<VECTOR<double,3>,1> const&,BOUNDARY_CONDITIONS_COLOR<VECTOR<double,3> >*,ARRAY<ARRAY<double> >&,int,double);
template void BASIS_INTEGRATION_UNIFORM_COLOR<VECTOR<double,3>,2>::Add_Surface_Block_Scalar<1>(SYSTEM_SURFACE_BLOCK_SCALAR_HELPER_COLOR<VECTOR<double,3> >&,
    BASIS_STENCIL_UNIFORM<VECTOR<double,3>,1> const&,BOUNDARY_CONDITIONS_SCALAR_COLOR<VECTOR<double,3> >*,ARRAY<ARRAY<double> >&,double);
template void BASIS_INTEGRATION_UNIFORM_COLOR<VECTOR<double,2>,2>::Add_Volume_Block<0,1>(SYSTEM_VOLUME_BLOCK_HELPER_COLOR<VECTOR<double,2> >&,
    BASIS_STENCIL_UNIFORM<VECTOR<double,2>,0> const&,BASIS_STENCIL_UNIFORM<VECTOR<double,2>,1> const&,ARRAY<double> const&);
template void BASIS_INTEGRATION_UNIFORM_COLOR<VECTOR<double,2>,2>::Add_Volume_Block<1,1>(SYSTEM_VOLUME_BLOCK_HELPER_COLOR<VECTOR<double,2> >&,
    BASIS_STENCIL_UNIFORM<VECTOR<double,2>,1> const&,BASIS_STENCIL_UNIFORM<VECTOR<double,2>,1> const&,ARRAY<double> const&);
template void BASIS_INTEGRATION_UNIFORM_COLOR<VECTOR<double,2>,2>::Add_Surface_Block<1>(SYSTEM_SURFACE_BLOCK_HELPER_COLOR<VECTOR<double,2> >&,
    BASIS_STENCIL_UNIFORM<VECTOR<double,2>,1> const&,BOUNDARY_CONDITIONS_COLOR<VECTOR<double,2> >*,ARRAY<ARRAY<double> >&,int,double);
template void BASIS_INTEGRATION_UNIFORM_COLOR<VECTOR<double,2>,2>::Add_Surface_Block_Scalar<1>(SYSTEM_SURFACE_BLOCK_SCALAR_HELPER_COLOR<VECTOR<double,2> >&,
    BASIS_STENCIL_UNIFORM<VECTOR<double,2>,1> const&,BOUNDARY_CONDITIONS_SCALAR_COLOR<VECTOR<double,2> >*,ARRAY<ARRAY<double> >&,double);
}
