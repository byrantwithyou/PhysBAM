//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Log/LOG.h>
#include <Core/Matrices/SPARSE_MATRIX_FLAT_MXN.h>
#include <Tools/Krylov_Solvers/CONJUGATE_GRADIENT.h>
#include <Tools/Krylov_Solvers/KRYLOV_VECTOR_WRAPPER.h>
#include <Tools/Krylov_Solvers/MATRIX_SYSTEM.h>
#include <Tools/Read_Write/OCTAVE_OUTPUT.h>
#include <Grid_Tools/Grids/FACE_RANGE_ITERATOR.h>
#include <Grid_PDE/Boundaries/BOUNDARY_MAC_GRID_PERIODIC.h>
#include <Geometry/Basic_Geometry/BASIC_SIMPLEX_POLICY.h>
#include <Geometry/Basic_Geometry/SEGMENT_2D.h>
#include <Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <Geometry/Grids_Uniform_Computations/MARCHING_CUBES.h>
#include <Geometry/Projection/CUT_CELL_PROJECTION.h>
namespace PhysBAM
{
//#####################################################################
// Function Cut_Cell_Projection
//#####################################################################
template<class TV> void CUT_CELL_PROJECTION<TV>::
Cut_Cell_Projection(const GRID<TV>& grid,int ghost,
    ARRAY<T,FACE_INDEX<TV::m> >& u,T density,T dt)
{
    typedef typename BASIC_SIMPLEX_POLICY<TV,TV::m>::SIMPLEX_FACE FACE;
    typedef KRYLOV_VECTOR_WRAPPER<T,ARRAY<T> > VEC;
    typedef MATRIX_SYSTEM<SPARSE_MATRIX_FLAT_MXN<T>,T,VEC> MAT;

    // TODO: how much mass does the surface velocity edge get?  mass -> 0?
    
    // TODO: fill in valid xfer data flags.

    static const int num_corners=1<<TV::m;
    static const int partial_object_flag=1<<29;
    static const int partial_surface_flag=1<<30;
    static const int index_mask=partial_object_flag-1;

    ARRAY<int,TV_INT> cell_index(grid.Domain_Indices(ghost),use_init,-9);
    VEC rhs,sol,null,precon_tmp;
    SPARSE_MATRIX_FLAT_MXN<T> A;
    MAT sys(A);
    
    // nonnegative = valid
    // -1 = air
    // -2 = object
    // -3 = domain boundary air
    // -4 = domain boundary object
    // -9 = uninitialized
    int next_cell_index=0;
    for(RANGE_ITERATOR<TV::m> it(grid.Domain_Indices());it.Valid();it.Next()){
        int& ci=cell_index(it.index);
        VECTOR<T,num_corners> phis;
        bool has_surface=false,has_object=false;
        if(surface_phi){
            MARCHING_CUBES<TV>::Compute_Phis_For_Cell(phis,*surface_phi,it.index);
            if(phis.Min()>=0){ci=-1;continue;} // all air
            has_surface=phis.Max()>0;}
        if(object_phi){
            MARCHING_CUBES<TV>::Compute_Phis_For_Cell(phis,*object_phi,it.index);
            if(phis.Min()>=0){ci=-1;continue;} // all air
            has_object=phis.Max()>0;}
        PHYSBAM_ASSERT(!has_object || !has_surface);

        ci=next_cell_index++;
        if(has_surface) ci|=partial_surface_flag;
        if(has_object) ci|=partial_object_flag;}
    rhs.v.Resize(next_cell_index,init_all,0);

    for(int s=0;s<2*TV::m;s++){
        if(bc_type(s)==0) continue;
        int value=-2-bc_type(s);
        for(RANGE_ITERATOR<TV::m> it(grid.Domain_Indices(),ghost,0,RI::ghost,s);it.Valid();it.Next())
            cell_index(it.index)=value;}

    VECTOR<bool,TV::m> is_periodic;
    for(int i=0;i<TV::m;i++) is_periodic(i)=!bc_type(2*i);
    if(bc_type.Contains(0))
        Fill_Ghost_Cells_Periodic(grid,cell_index,cell_index,is_periodic,ghost);
    
    ARRAY<VECTOR<TV,TV::m> > surface;
    VECTOR<ARRAY<VECTOR<TV,TV::m> >,2*TV::m> boundary;
    VECTOR<VECTOR<ARRAY<VECTOR<TV,TV::m> >*,2>,2*TV::m> pboundary;
    for(int s=0;s<2*TV::m;s++) pboundary(s)(1)=&boundary(s);
    
    // Assumes square cells
    T one_over_dx=grid.one_over_dX(0);
    TV base_entry=grid.one_over_dX*grid.one_over_dX/density;
    // times 2 to compensate for pressure at face
    T base_d_entry=one_over_dx*one_over_dx/density*2;
    TV sc=grid.one_over_dX/density;
    bool has_nullspace=true;

    A.Reset(next_cell_index);
    for(RANGE_ITERATOR<TV::m> it(grid.Domain_Indices());it.Valid();it.Next())
    {
        int this_index=cell_index(it.index);
        if(this_index<0) continue;
        bool has_surface=this_index&partial_surface_flag;
        bool has_object=this_index&partial_object_flag;
        this_index&=index_mask;
        
        VECTOR<T,TV::m*2> face_fraction;
        TV surface_area_weighted_normal;
        TV centroid;
        T total_area=0;
        T surface_fraction=0;
        
        T diag_entry=0;
        if(has_surface || has_object){
            VECTOR<T,num_corners> phis;
            if(has_surface)
                MARCHING_CUBES<TV>::Compute_Phis_For_Cell(phis,*surface_phi,it.index);
            if(has_object)
                MARCHING_CUBES<TV>::Compute_Phis_For_Cell(phis,*object_phi,it.index);
            surface.Remove_All();
            for(int s=0;s<2*TV::m;s++) boundary(s).Remove_All();
            MARCHING_CUBES<TV>::Get_Elements_For_Cell(surface,pboundary,phis);
            for(int s=0;s<2*TV::m;s++)
                for(const auto& t:boundary(s))
                    face_fraction(s)+=FACE::Size(t);
            for(const auto& t:surface)
            {
                TV an=FACE::Area_Weighted_Normal(t);
                T area=an.Magnitude();
                surface_area_weighted_normal+=an;
                centroid+=t.Sum()/t.m*area;
                total_area+=area;
            }
            centroid/=total_area;
            centroid=grid.Cell_Domain(it.index).Point_From_Normalized_Coordinates(centroid);
            PHYSBAM_ASSERT(grid.Cell_Domain(it.index).Lazy_Inside(centroid));
            if(has_surface)
            {
                if(bc_p){
                    T p=bc_p(centroid)*dt;
                    rhs.v(this_index)+=base_d_entry*p*surface_area_weighted_normal.Magnitude();}
                // TODO: need rhs velocity at surface edge
                has_nullspace=false;
                diag_entry+=base_d_entry*surface_area_weighted_normal.Magnitude();
            }
            else
            {
                if(bc_v){
                    TV v=bc_v(centroid);
                    rhs.v(this_index)-=one_over_dx*surface_area_weighted_normal.Dot(v);}
            }
        }
        else face_fraction.Fill(1);
        
        for(int s=0;s<2*TV::m;s++)
        {
            if(!face_fraction(s)) continue;
            FACE_INDEX<TV::m> face(s/2,it.index);
            face.index(s/2)+=s%2;
            T sign=s%2?1:-1;
            TV_INT cell(it.index);
            cell(s/2)+=sign;
            int ci=cell_index(cell);
            PHYSBAM_ASSERT(ci>=0 || ci==-4 || ci==-3 || !face_fraction(s));
            T entry=base_entry(s/2)*face_fraction(s);
            T div_entry=sign*grid.one_over_dX(s/2)*face_fraction(s);
            if(ci==-4)
            {
                T v=bc_v?bc_v(grid.Face(face))(s/2):0;
                if(!has_surface) rhs.v(this_index)-=div_entry*v;
                u(face)=v;
                continue;
            }
            if(!has_surface) rhs.v(this_index)-=div_entry*u(face);
            if(ci==-3)
            {
                entry*=2; // compensate for pressure at face
                if(bc_p){
                    T p=bc_p(grid.Face(face))*dt;
                    rhs.v(this_index)+=entry*p;
                    u(face)-=sign*sc(s/2)*p;}
                has_nullspace=false;
            }
            else
            {
                A.Append_Entry_To_Current_Row(ci&index_mask,-entry);
            }
            diag_entry+=entry;
        }
        A.Append_Entry_To_Current_Row(this_index,diag_entry);
        A.Finish_Row();
    }
    A.Sort_Entries();
    sol.v.Resize(rhs.v.m,init_all,0);
    if(has_nullspace){
        null.v.Resize(rhs.v.m,init_all,sqrt((T)1/rhs.v.m));
        sys.nullspace_vectors.Append(&null);}

    sys.use_preconditioner=use_preconditioner;
    sys.preconditioner_commutes_with_projection=(!has_nullspace||!use_preconditioner);
    if(use_preconditioner){
        A.Construct_Incomplete_Cholesky_Factorization();
        sys.P=A.C;
        precon_tmp.v.Resize(next_cell_index,init_all,0);
        sys.temp_vector=&precon_tmp;}
    
    static int solve_id=-1;
    solve_id++;
    if(test_system) sys.Test_System(sol);
    if(print_matrix){
        LOG::cout<<"solve id: "<<solve_id<<std::endl;
        OCTAVE_OUTPUT<T>(LOG::sprintf("M-%i.txt",solve_id).c_str()).Write("M",sys,rhs);
        OCTAVE_OUTPUT<T>(LOG::sprintf("C-%i.txt",solve_id).c_str()).Write_Preconditioner("C",sys,rhs);
        OCTAVE_OUTPUT<T>(LOG::sprintf("P-%i.txt",solve_id).c_str()).Write_Projection("P",sys,rhs);
        OCTAVE_OUTPUT<T>(LOG::sprintf("b-%i.txt",solve_id).c_str()).Write("b",rhs);}

    ARRAY<KRYLOV_VECTOR_BASE<T>*> av;
    CONJUGATE_GRADIENT<T> cg;
    cg.finish_before_indefiniteness=true;
    cg.relative_tolerance=true;
    bool converged=cg.Solve(sys,sol,rhs,av,solver_tolerance,0,solver_iterations);
    if(!converged) LOG::printf("SOLVER DID NOT CONVERGE.\n");

    if(print_matrix)
        OCTAVE_OUTPUT<T>(LOG::sprintf("x-%i.txt",solve_id).c_str()).Write("x",sol);

    if(print_residual){
        sys.Multiply(sol,*av(0));
        *av(0)-=rhs;
        T r=sys.Convergence_Norm(*av(0));
        LOG::cout<<"residual: "<<r<<std::endl;}

    if(valid_u) valid_u->Resize(grid.Domain_Indices(ghost),false,false);

    for(FACE_RANGE_ITERATOR<TV::m> it(grid.Domain_Indices(ghost),RF::skip_outer);it.Valid();it.Next()){
        int ca=cell_index(it.face.First_Cell_Index());
        int cb=cell_index(it.face.Second_Cell_Index());
        if((ca<0 && cb<0) || ca==-1 || ca==-2 || cb==-1 || cb==-2){
            if(valid_u) (*valid_u)(it.face)=false;
            continue;}
        if(valid_u) (*valid_u)(it.face)=true;
        if(ca==-4 || cb==-4) continue;
        T pa=ca<0?0:sol.v(ca&index_mask),pb=cb<0?0:sol.v(cb&index_mask);
        u(it.face)-=(pb-pa)*sc(it.face.axis);}

    if(valid_p || pressure){
        if(valid_p) valid_p->Resize(grid.Domain_Indices(ghost),no_init);
        if(pressure) pressure->Resize(grid.Domain_Indices(ghost),no_init);
        for(RANGE_ITERATOR<TV::m> it(valid_p->domain);it.Valid();it.Next()){
            int c=cell_index(it.index);
            if(c<0) continue;
            if(valid_p) (*valid_p)(it.index)=true;
            if(pressure) (*pressure)(it.index)=sol.v(c&index_mask)/dt;}}

    if(bc_type.Contains(0)){
        Fill_Ghost_Faces_Periodic(grid,u,u,is_periodic,ghost);
        if(valid_u) Fill_Ghost_Faces_Periodic(grid,*valid_u,*valid_u,is_periodic,ghost);}
}
//#####################################################################
template class CUT_CELL_PROJECTION<VECTOR<float,2> >;
template class CUT_CELL_PROJECTION<VECTOR<float,3> >;
template class CUT_CELL_PROJECTION<VECTOR<double,2> >;
template class CUT_CELL_PROJECTION<VECTOR<double,3> >;
}
