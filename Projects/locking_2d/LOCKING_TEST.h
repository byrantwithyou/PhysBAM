//#####################################################################
// Copyright 2012, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef __LOCKING_TEST__
#define __LOCKING_TEST__

#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Matrices/MATRIX_MXN.h>

#define U_GRID(X) (X)%n
#define V_GRID(X) (X)%n+n
#define P_GRID(X) (X)%n+2*n

namespace PhysBAM{

template<class T>
class LOCKING_TEST:public NONCOPYABLE
{
private:
    
    bool affine_velicities; // otherwise bilinear
    int n; // cells
    T mu; // viscosity

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

    LOCKING_TEST():affine_velicities(false),n(2),mu(3){}
    ~LOCKING_TEST(){};

    void Parse_Arguments(int argc,char *argv[])
    {
        PARSE_ARGS parse_args;
        parse_args.Add_Option_Argument("-affine","use affine velocities");
        parse_args.Add_Integer_Argument("-cells",2,"cells","number of grid cells");
        parse_args.Add_Double_Argument("-viscosity",3,"viscosity");
        parse_args.Parse(argc,argv);
        affine_velicities=parse_args.Is_Value_Set("-affine");
        n=parse_args.Get_Integer_Value("-cells");
        mu=parse_args.Get_Double_Value("-viscosity");
    }

    void Compute_System_Matrix()
    {
        Generate_Element_Stiffness_Matrices();
        A.Resize(n*3,n*3);
        for(int x=0; x<n; x++) for(int y=0; y<n; y++)
        {
            for(int ix=0; ix<2; ix++) for(int iy=0; iy<2; iy++) 
            for(int jx=0; jx<2; jx++) for(int jy=0; jy<2; jy++)
            {
                A(U_GRID(x+ix),U_GRID(y+iy))+=(2*NxNx[ix][iy][jx][jy]+NyNy[ix][iy][jx][jy])*mu;
                A(V_GRID(x+ix),V_GRID(y+iy))+=(2*NyNy[ix][iy][jx][jy]+NxNx[ix][iy][jx][jy])*mu;
                A(V_GRID(x+ix),U_GRID(y+iy))+=(NxNy[ix][iy][jx][jy])*mu;
                A(U_GRID(y+iy),V_GRID(x+ix))+=(NxNy[ix][iy][jx][jy])*mu;
            }
            for(int ix=0; ix<2; ix++) for(int iy=0; iy<2; iy++)
            {
                
            }
        }
    }
};
}
#endif
