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
#include <PhysBAM_Tools/Read_Write/Octave/OCTAVE_OUTPUT.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Basic_Geometry/SPHERE.h>
#include <PhysBAM_Geometry/Finite_Elements/BASIS_INTEGRATION_UNIFORM.h>
#include <PhysBAM_Geometry/Finite_Elements/BASIS_STENCIL_UNIFORM.h>
#include <PhysBAM_Geometry/Finite_Elements/CELL_MAPPING.h>
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

    GRID<TV> grid(TV_INT()+8,RANGE<TV>(TV(),TV()+1),true);
    GRID<TV> coarse_grid(grid.counts/2,grid.domain,true);
    ARRAY<T,TV_INT> phi(coarse_grid.Node_Indices());
    SPHERE<TV> sphere(TV()+(T).5,(T).5);
    for(UNIFORM_GRID_ITERATOR_NODE<TV> it(coarse_grid);it.Valid();it.Next())
        phi(it.index)=sphere.Signed_Distance(coarse_grid.Node(it.index));

    BASIS_STENCIL_UNIFORM<TV> p_stencil,u_stencil,v_stencil;
    p_stencil.Set_Center();
    p_stencil.Set_Constant_Stencil();
    p_stencil.Dice_Stencil();
    u_stencil.Set_Face(0);
    u_stencil.Set_Multilinear_Stencil();
    u_stencil.Dice_Stencil();
    v_stencil.Set_Face(1);
    v_stencil.Set_Multilinear_Stencil();
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

    CELL_MAPPING<TV> index_map_p(grid),index_map_u(grid),index_map_v(grid);
    index_map_p.periodic.Fill(true);
    index_map_u.periodic.Fill(true);
    index_map_v.periodic.Fill(true);

    SYSTEM_MATRIX_HELPER<T> helper;
    RANGE<TV_INT> boundary_conditions;
    boundary_conditions.min_corner.Fill(BASIS_INTEGRATION_UNIFORM<TV>::periodic);
    boundary_conditions.max_corner.Fill(BASIS_INTEGRATION_UNIFORM<TV>::periodic);

    LOG::cout<<"u ";u_stencil.Print();
    LOG::cout<<"v ";v_stencil.Print();
    LOG::cout<<"p ";p_stencil.Print();
    LOG::cout<<"udx ";udx_stencil.Print();
    LOG::cout<<"vdx ";vdx_stencil.Print();
    LOG::cout<<"udy ";udy_stencil.Print();
    LOG::cout<<"vdy ";vdy_stencil.Print();

    LOG::cout<<"udx udx"<<std::endl;
    BASIS_INTEGRATION_UNIFORM<TV>(boundary_conditions,grid,grid,udx_stencil,udx_stencil,index_map_u,index_map_u,phi).Compute_Matrix(helper);
    helper.Scale(2*mu);
    LOG::cout<<"udy udy"<<std::endl;
    BASIS_INTEGRATION_UNIFORM<TV>(boundary_conditions,grid,grid,udy_stencil,udy_stencil,index_map_u,index_map_u,phi).Compute_Matrix(helper);
    helper.Scale(mu);
    LOG::cout<<"udy vdx"<<std::endl;
    BASIS_INTEGRATION_UNIFORM<TV>(boundary_conditions,grid,grid,udy_stencil,vdx_stencil,index_map_u,index_map_v,phi).Compute_Matrix(helper);
    helper.Scale(mu);
    LOG::cout<<"vdx udy"<<std::endl;
    BASIS_INTEGRATION_UNIFORM<TV>(boundary_conditions,grid,grid,vdx_stencil,udy_stencil,index_map_v,index_map_u,phi).Compute_Matrix(helper);
    helper.Scale(mu);
    LOG::cout<<"vdx vdx"<<std::endl;
    BASIS_INTEGRATION_UNIFORM<TV>(boundary_conditions,grid,grid,vdx_stencil,vdx_stencil,index_map_v,index_map_v,phi).Compute_Matrix(helper);
    helper.Scale(2*mu);
    LOG::cout<<"vdy vdy"<<std::endl;
    BASIS_INTEGRATION_UNIFORM<TV>(boundary_conditions,grid,grid,vdy_stencil,vdy_stencil,index_map_v,index_map_v,phi).Compute_Matrix(helper);
    helper.Scale(mu);
    LOG::cout<<"udx p"<<std::endl;
    BASIS_INTEGRATION_UNIFORM<TV>(boundary_conditions,grid,grid,udx_stencil,p_stencil,index_map_u,index_map_p,phi).Compute_Matrix(helper);
    LOG::cout<<"vdy p"<<std::endl;
    BASIS_INTEGRATION_UNIFORM<TV>(boundary_conditions,grid,grid,vdy_stencil,p_stencil,index_map_v,index_map_p,phi).Compute_Matrix(helper);
    LOG::cout<<"p udx"<<std::endl;
    BASIS_INTEGRATION_UNIFORM<TV>(boundary_conditions,grid,grid,p_stencil,udx_stencil,index_map_p,index_map_u,phi).Compute_Matrix(helper);
    LOG::cout<<"p vdy"<<std::endl;
    BASIS_INTEGRATION_UNIFORM<TV>(boundary_conditions,grid,grid,p_stencil,vdy_stencil,index_map_p,index_map_v,phi).Compute_Matrix(helper);
    SPARSE_MATRIX_FLAT_MXN<T> matrix;
    helper.Set_Matrix(3*grid.counts.Product(),3*grid.counts.Product(),matrix,1e-14);

    OCTAVE_OUTPUT<T>("M.txt").Write("M",matrix);
};

int main(int argc,char* argv[])
{
    Integration_Test(argc,argv);

    return 0;
}
