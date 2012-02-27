//#####################################################################
// Copyright 2012, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef __LOCKING_TEST__
#define __LOCKING_TEST__

#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Matrices/MATRIX_MXN.h>

namespace PhysBAM{

template<class T>
class LOCKING_TEST:public NONCOPYABLE
{
private:
    
    bool affine_velicities; // otherwise bilinear
    long cells;

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
        }
        else
        {
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

    LOCKING_TEST():affine_velicities(false),cells(2){}
    ~LOCKING_TEST(){};

    void Parse_Arguments(int argc,char *argv[])
    {
        PARSE_ARGS parse_args;
        parse_args.Add_Option_Argument("-affine","use affine velocities");
        parse_args.Add_Integer_Argument("-cells",2,"cells","number of grid cells");
        parse_args.Parse(argc,argv);
        affine_velicities=parse_args.Is_Value_Set("-affine");
        cells=parse_args.Get_Integer_Value("-cells");
    }

    void Compute_System_Matrix()
    {
        Generate_Element_Stiffness_Matrices();
        A.Resize(cells*3,cells*3);
        
    }
};
}
#endif
