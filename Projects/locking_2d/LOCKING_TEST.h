//#####################################################################
// Copyright 2012, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef __LOCKING_TEST__
#define __LOCKING_TEST__

#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Matrices/MATRIX_MXN.h>
#include <PhysBAM_Tools/Matrices/SPARSE_MATRIX_FLAT_MXN.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Read_Write/Octave/OCTAVE_OUTPUT.h>
#include <PhysBAM_Geometry/Finite_Elements/BASIS_INTEGRATION_UNIFORM.h>
#include <PhysBAM_Geometry/Finite_Elements/BASIS_STENCIL_UNIFORM.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/SYSTEM_MATRIX_HELPER.h>

#define U_GRID(X,Y) ((X)%n)*n+((Y)%n)
#define V_GRID(X,Y) ((X)%n)*n+((Y)%n)+N
#define P_GRID(X,Y) ((X)%n)*n+((Y)%n)+2*N

namespace PhysBAM{

template<class T>
class LOCKING_TEST:public NONCOPYABLE
{
private:
    
    bool affine_velicities; // otherwise bilinear
    int n; // cells per direction
    int N; // cells
    T mu; // viscosity
    bool use_framework;

    MATRIX_MXN<T> A;

    T NxNx[2][2][2][2];
    T NyNy[2][2][2][2];
    T NxNy[2][2][2][2];
    T NxPhi[2][2];
    T NyPhi[2][2];
    
    void Generate_Element_Stiffness_Matrices()
    {
        if(affine_velicities)
        {
            for(int ix=0; ix<2; ix++) for(int iy=0; iy<2; iy++) 
            for(int jx=0; jx<2; jx++) for(int jy=0; jy<2; jy++)
            { 
                NxNx[ix][iy][jx][jy]=((ix==jx)?(0.25):(-0.25));
                NyNy[ix][iy][jx][jy]=((iy==jy)?(0.25):(-0.25));
                NxNy[ix][iy][jx][jy]=((ix==jy)?(0.25):(-0.25));
            }
            for(int ix=0; ix<2; ix++) for(int iy=0; iy<2; iy++)
            {
                NxPhi[ix][iy]=((ix>0)?(0.5):(-0.5));
                NyPhi[ix][iy]=((iy>0)?(0.5):(-0.5));
            }
        }else{
            for(int ix=0; ix<2; ix++) for(int iy=0; iy<2; iy++) 
            for(int jx=0; jx<2; jx++) for(int jy=0; jy<2; jy++)
            { 
                NxNx[ix][iy][jx][jy]=((ix==jx)?(1):(-1))*((iy==jy)?(1):(0.5))/3.0;
                NyNy[ix][iy][jx][jy]=((iy==jy)?(1):(-1))*((ix==jx)?(1):(0.5))/3.0;
                NxNy[ix][iy][jx][jy]=((ix==jy)?(1):(-1))*0.25;
            }
            for(int ix=0; ix<2; ix++) for(int iy=0; iy<2; iy++)
            {
                NxPhi[ix][iy]=((ix>0)?(0.5):(-0.5));
                NyPhi[ix][iy]=((iy>0)?(0.5):(-0.5));
            }
        }
    }

public:

    LOCKING_TEST():affine_velicities(false),n(2),N(n*n),mu(3),use_framework(false) {}
    ~LOCKING_TEST(){};

    void Parse_Arguments(int argc,char *argv[])
    {
        PARSE_ARGS parse_args;
        parse_args.Add_Option_Argument("-affine","use affine velocities");
        parse_args.Add_Option_Argument("-use_framework","use integration framework");
        parse_args.Add_Integer_Argument("-cells",2,"cells","number of grid cells");
        parse_args.Add_Double_Argument("-viscosity",3,"viscosity");
        parse_args.Parse(argc,argv);
        affine_velicities=parse_args.Is_Value_Set("-affine");
        n=parse_args.Get_Integer_Value("-cells"); N=n*n;
        mu=parse_args.Get_Double_Value("-viscosity");
        use_framework=parse_args.Is_Value_Set("-use_framework");

        if(!use_framework) return;
        typedef VECTOR<int,2> TV_INT;
        typedef VECTOR<T,2> TV;
        GRID<TV> grid(TV_INT(n,n),RANGE<TV>::Unit_Box(),true);

        BASIS_STENCIL_UNIFORM<TV> p_stencil,u_stencil,v_stencil;
        p_stencil.Set_Constant_Stencil();
        p_stencil.Set_Center();
        u_stencil.Set_Multilinear_Stencil();
        u_stencil.Set_Face(0);
        v_stencil.Set_Multilinear_Stencil();
        v_stencil.Set_Face(1);
        BASIS_STENCIL_UNIFORM<TV> udx_stencil(u_stencil),udy_stencil(u_stencil),vdx_stencil(v_stencil),vdy_stencil(v_stencil);
        udx_stencil.Differentiate(0);
        udy_stencil.Differentiate(1);
        vdx_stencil.Differentiate(0);
        vdy_stencil.Differentiate(1);

        ARRAY<int,TV_INT> index_map_p(grid.Domain_Indices());
        int k=0;
        for(UNIFORM_GRID_ITERATOR_CELL<TV> it(grid);it.Valid();it.Next())
            index_map_p(it.index)=k++;
        ARRAY<int,TV_INT> index_map_u(index_map_p),index_map_v(index_map_p);
        index_map_u+=grid.counts.Product();
        index_map_v+=2*grid.counts.Product();

        SYSTEM_MATRIX_HELPER<T> helper;
        BASIS_INTEGRATION_UNIFORM<TV> biu(grid);
        biu.boundary_conditions.min_corner.Fill(BASIS_INTEGRATION_UNIFORM<TV>::periodic);
        biu.boundary_conditions.max_corner.Fill(BASIS_INTEGRATION_UNIFORM<TV>::periodic);

        biu.Compute_Matrix(helper, udx_stencil, udx_stencil, index_map_u, index_map_u);
        helper.Scale(2*mu);
        biu.Compute_Matrix(helper, udy_stencil, udy_stencil, index_map_u, index_map_u);
        helper.Scale(mu);
        biu.Compute_Matrix(helper, udy_stencil, vdx_stencil, index_map_u, index_map_v);
        helper.Scale(mu);
        biu.Compute_Matrix(helper, vdx_stencil, udy_stencil, index_map_v, index_map_u);
        helper.Scale(mu);
        biu.Compute_Matrix(helper, vdx_stencil, vdx_stencil, index_map_v, index_map_v);
        helper.Scale(2*mu);
        biu.Compute_Matrix(helper, vdy_stencil, vdy_stencil, index_map_v, index_map_v);
        helper.Scale(mu);
        biu.Compute_Matrix(helper, udx_stencil, p_stencil, index_map_u, index_map_p);
        helper.Scale(mu);
        biu.Compute_Matrix(helper, vdy_stencil, p_stencil, index_map_v, index_map_p);
        helper.Scale(mu);
        SPARSE_MATRIX_FLAT_MXN<T> matrix;
        helper.Set_Matrix(3*grid.counts.Product(),3*grid.counts.Product(),matrix);

        OCTAVE_OUTPUT<T>("M.txt").Write("M",matrix);
    }

    void Compute_System_Matrix()
    {
        Generate_Element_Stiffness_Matrices();
        A.Resize(N*3,N*3);
        for(int x=0; x<n; x++) for(int y=0; y<n; y++)
        {
            for(int ix=0; ix<2; ix++) for(int iy=0; iy<2; iy++) 
            for(int jx=0; jx<2; jx++) for(int jy=0; jy<2; jy++)
            {
                A(U_GRID(x+ix,y+iy),U_GRID(x+jx,y+jy))+=(2*NxNx[ix][iy][jx][jy]+NyNy[ix][iy][jx][jy])*mu;
                A(V_GRID(x+ix,y+iy),V_GRID(x+jx,y+jy))+=(2*NyNy[ix][iy][jx][jy]+NxNx[ix][iy][jx][jy])*mu;
                A(V_GRID(x+ix,y+iy),U_GRID(x+jx,y+jy))+=(NxNy[ix][iy][jx][jy])*mu;
                A(U_GRID(x+jx,y+jy),V_GRID(x+ix,y+iy))+=(NxNy[ix][iy][jx][jy])*mu;
            }
            for(int ix=0; ix<2; ix++) for(int iy=0; iy<2; iy++)
            {
                A(U_GRID(x+ix,y+iy),P_GRID(x,y))-=NxPhi[ix][iy];
                A(V_GRID(x+ix,y+iy),P_GRID(x,y))-=NyPhi[ix][iy];
                A(P_GRID(x,y),U_GRID(x+ix,y+iy))-=NxPhi[ix][iy];
                A(P_GRID(x,y),V_GRID(x+ix,y+iy))-=NyPhi[ix][iy];
            }
        }
    }

    void Print_System_Matrix()
    {
        for(int x=0; x<3*N; x++)
        {
            for(int y=0; y<3*N; y++)
                LOG::cout<<A(x,y)<<"\t";
            LOG::cout<<std::endl;
        }                            
    }
};
}
#endif
