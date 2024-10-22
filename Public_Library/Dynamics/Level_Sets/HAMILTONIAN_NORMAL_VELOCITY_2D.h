//#####################################################################
// Copyright 2002, Ronald Fedkiw, Frederic Gibou.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class HAMILTONIAN_NORMAL_VELOCITY_2D 
//##################################################################### 
#ifndef __HAMILTONIAN_NORMAL_VELOCITY_2D__
#define __HAMILTONIAN_NORMAL_VELOCITY_2D__   

#include <Core/Arrays/ARRAYS_FORWARD.h>
#include <Dynamics/Level_Sets/HAMILTONIAN_2D.h>
namespace PhysBAM{

template<class T>
class HAMILTONIAN_NORMAL_VELOCITY_2D:public HAMILTONIAN_2D<T>
{
    typedef VECTOR<T,2> TV;
public:
    ARRAY<T,VECTOR<int,2> >& speed;

    HAMILTONIAN_NORMAL_VELOCITY_2D(GRID<TV>& grid_input,ARRAY<T,VECTOR<int,2> >& speed_input)
        :HAMILTONIAN_2D<T>(grid_input),speed(speed_input)
    {}

//#####################################################################
    T H(const T phi_x,const T phi_y,const int i=0,const int j=0,const T t=0) override;
    T Maxabs_H1(const T phi_x_1,const T phi_x_2,const T phi_y_1,const T phi_y_2,const int i=0,const int j=0,const T t=0) override;
    T Maxabs_H2(const T phi_x_1,const T phi_x_2,const T phi_y_1,const T phi_y_2,const int i=0,const int j=0,const T t=0) override;
//#####################################################################
};   
}
#endif
