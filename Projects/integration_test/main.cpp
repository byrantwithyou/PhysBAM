//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Matrices/SPARSE_MATRIX_FLAT_MXN.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Read_Write/OCTAVE_OUTPUT.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Basic_Geometry/SPHERE.h>
#include <PhysBAM_Geometry/Finite_Elements/BASIS_INTEGRATION_BOUNDARY_UNIFORM.h>
#include <PhysBAM_Geometry/Finite_Elements/BASIS_INTEGRATION_UNIFORM.h>
#include <PhysBAM_Geometry/Finite_Elements/BASIS_STENCIL_BOUNDARY_UNIFORM.h>
#include <PhysBAM_Geometry/Finite_Elements/BASIS_STENCIL_UNIFORM.h>
#include <PhysBAM_Geometry/Finite_Elements/CELL_MAPPING.h>
#include <PhysBAM_Geometry/Grids_Uniform_Computations/MARCHING_CUBES.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/SYSTEM_MATRIX_HELPER.h>

using namespace PhysBAM;

typedef double T;
typedef VECTOR<T,2> TV;
typedef VECTOR<int,2> TV_INT;
typedef float RW;

// * Handle interface
// * Handle other boundary conditions
// * Handle 

void Integration_Test(int argc,char* argv[])
{
    PARSE_ARGS parse_args;
    parse_args.Add_Double_Argument("-viscosity",1,"viscosity");
    parse_args.Add_Double_Argument("-m",1,"meter scale");
    parse_args.Add_Double_Argument("-s",1,"second scale");
    parse_args.Add_Double_Argument("-kg",1,"kilogram scale");
    parse_args.Parse(argc,argv);
    T m=parse_args.Get_Double_Value("-m");
    T s=parse_args.Get_Double_Value("-s");
    T kg=parse_args.Get_Double_Value("-kg");
    T mu=parse_args.Get_Double_Value("-viscosity")*kg/s;
    (void)m;

    GRID<TV> grid(TV_INT()+4,RANGE<TV>(TV(),TV()+1)*m,true);
    GRID<TV> coarse_grid(grid.counts,grid.domain,true);
    ARRAY<T,TV_INT> phi(coarse_grid.Node_Indices());
    SPHERE<TV> sphere(TV()+(T).5*m,(T).5*m);
    for(UNIFORM_GRID_ITERATOR_NODE<TV> it(coarse_grid);it.Valid();it.Next())
        phi(it.index)=sphere.Signed_Distance(coarse_grid.Node(it.index));

    BASIS_STENCIL_UNIFORM<TV> p_stencil,u_stencil,v_stencil;
    p_stencil.Set_Center();
    p_stencil.Set_Constant_Stencil();
    p_stencil.Scale_Axes(grid.one_over_dX);
    p_stencil.Dice_Stencil();
    u_stencil.Set_Face(0);
    u_stencil.Set_Multilinear_Stencil();
    u_stencil.Scale_Axes(grid.one_over_dX);
    u_stencil.Dice_Stencil();
    v_stencil.Set_Face(1);
    v_stencil.Set_Multilinear_Stencil();
    v_stencil.Scale_Axes(grid.one_over_dX);
    v_stencil.Dice_Stencil();
    BASIS_STENCIL_UNIFORM<TV> udx_stencil(u_stencil),udy_stencil(u_stencil),vdx_stencil(v_stencil),vdy_stencil(v_stencil);
    udx_stencil.Differentiate(0);
    udx_stencil.Dice_Stencil();
    udy_stencil.Differentiate(1);
    udy_stencil.Dice_Stencil();
    vdx_stencil.Differentiate(0);
    vdx_stencil.Dice_Stencil();
    vdy_stencil.Differentiate(1);
    vdy_stencil.Dice_Stencil();

    SEGMENTED_CURVE_2D<T> curve;
    MARCHING_CUBES<TV>::Create_Surface(curve,coarse_grid,phi);
    BASIS_STENCIL_BOUNDARY_UNIFORM<TV> q_stencil(curve);

    CELL_MAPPING<TV> index_map_p(grid),index_map_u(grid),index_map_v(grid);
    index_map_p.periodic.Fill(true);
    index_map_u.periodic.Fill(true);
    index_map_v.periodic.Fill(true);

    SYSTEM_MATRIX_HELPER<T> helper;
    RANGE<TV_INT> boundary_conditions;
    boundary_conditions.min_corner.Fill(BASIS_INTEGRATION_UNIFORM<TV>::periodic);
    boundary_conditions.max_corner.Fill(BASIS_INTEGRATION_UNIFORM<TV>::periodic);

    BASIS_INTEGRATION_UNIFORM<TV>(boundary_conditions,grid,coarse_grid,udx_stencil,udx_stencil,index_map_u,index_map_u,phi).Compute_Matrix(helper);
    helper.Scale(2*mu);
    BASIS_INTEGRATION_UNIFORM<TV>(boundary_conditions,grid,coarse_grid,udy_stencil,udy_stencil,index_map_u,index_map_u,phi).Compute_Matrix(helper);
    helper.Scale(mu);
    helper.start=0;
    INTERVAL<int> block_uu=helper.Get_Block();

    BASIS_INTEGRATION_UNIFORM<TV>(boundary_conditions,grid,coarse_grid,vdx_stencil,vdx_stencil,index_map_v,index_map_v,phi).Compute_Matrix(helper);
    helper.Scale(2*mu);
    BASIS_INTEGRATION_UNIFORM<TV>(boundary_conditions,grid,coarse_grid,vdy_stencil,vdy_stencil,index_map_v,index_map_v,phi).Compute_Matrix(helper);
    helper.Scale(mu);
    helper.start=block_uu.max_corner;
    INTERVAL<int> block_vv=helper.Get_Block();

    BASIS_INTEGRATION_UNIFORM<TV>(boundary_conditions,grid,coarse_grid,udy_stencil,vdx_stencil,index_map_u,index_map_v,phi).Compute_Matrix(helper);
    helper.Scale(mu);
    INTERVAL<int> block_uv=helper.Get_Block();

    BASIS_INTEGRATION_UNIFORM<TV>(boundary_conditions,grid,coarse_grid,udx_stencil,p_stencil,index_map_u,index_map_p,phi).Compute_Matrix(helper);
    INTERVAL<int> block_up=helper.Get_Block();

    BASIS_INTEGRATION_UNIFORM<TV>(boundary_conditions,grid,coarse_grid,vdy_stencil,p_stencil,index_map_v,index_map_p,phi).Compute_Matrix(helper);
    INTERVAL<int> block_vp=helper.Get_Block();

    BASIS_INTEGRATION_BOUNDARY_UNIFORM<TV>(grid,u_stencil,q_stencil,index_map_u).Compute_Matrix(helper);
    INTERVAL<int> block_uq=helper.Get_Block();

    BASIS_INTEGRATION_BOUNDARY_UNIFORM<TV>(grid,v_stencil,q_stencil,index_map_v).Compute_Matrix(helper);
    INTERVAL<int> block_vq=helper.Get_Block();

    int start_v=index_map_u.next_index;
    int start_p=index_map_v.next_index+start_v;
    int start_uq=index_map_p.next_index+start_p;
    int start_vq=curve.mesh.elements.m+start_uq;
    int total=curve.mesh.elements.m+start_vq;

    printf("u %i v %i p %i uq %i vq %i end %i\n", 0, start_v, start_p, start_uq, start_vq, total);

    helper.Shift(start_v,start_v,block_vv);
    helper.Shift(0,start_v,block_uv);
    helper.Shift(0,start_p,block_up);
    helper.Shift(start_v,start_p,block_vp);
    helper.Shift(0,start_uq,block_uq);
    helper.Shift(start_v,start_vq,block_vq);
    helper.Add_Transpose(INTERVAL<int>(block_uv.min_corner,block_vq.max_corner));

    SPARSE_MATRIX_FLAT_MXN<T> matrix;
    helper.Set_Matrix(total,total,matrix,1e-14);

    OCTAVE_OUTPUT<T>("M.txt").Write("M",matrix);
};

int main(int argc,char* argv[])
{
    Integration_Test(argc,argv);

    return 0;
}
