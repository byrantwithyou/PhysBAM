//#####################################################################
// Copyright 2002, Ronald Fedkiw.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class HAMILTONIAN_NORMAL_ADVECTION_3D 
//##################################################################### 
#ifndef __HAMILTONIAN_NORMAL_ADVECTION_3D__
#define __HAMILTONIAN_NORMAL_ADVECTION_3D__   

#include <Dynamics/Level_Sets/HAMILTONIAN_3D.h>
namespace PhysBAM{

template<class T>
class HAMILTONIAN_NORMAL_ADVECTION_3D:public HAMILTONIAN_3D<T>
{
    typedef VECTOR<T,3> TV;
public:
    T speed;

    HAMILTONIAN_NORMAL_ADVECTION_3D(GRID<TV>& grid_input)
        :HAMILTONIAN_3D<T>(grid_input)
    {
        Set_Speed(1);
    }

    void Set_Speed(const T speed_input)
    {speed=speed_input;}

//#####################################################################
    T H(const T phi_x,const T phi_y,const T phi_z,const int i=0,const int j=0,const int k=0,const T t=0) override;
    T Maxabs_H1(const T phi_x_1,const T phi_x_2,const T phi_y_1,const T phi_y_2,const T phi_z_1,const T phi_z_2,const int i=0,const int j=0,const int k=0,const T t=0) override;
    T Maxabs_H2(const T phi_x_1,const T phi_x_2,const T phi_y_1,const T phi_y_2,const T phi_z_1,const T phi_z_2,const int i=0,const int j=0,const int k=0,const T t=0) override;
    T Maxabs_H3(const T phi_x_1,const T phi_x_2,const T phi_y_1,const T phi_y_2,const T phi_z_1,const T phi_z_2,const int i=0,const int j=0,const int k=0,const T t=0) override;
//#####################################################################
};   
}
#endif

