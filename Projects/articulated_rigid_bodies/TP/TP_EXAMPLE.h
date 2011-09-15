//#####################################################################
// Copyright 2004, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TP_EXAMPLE
//##################################################################### 
#ifndef __TP_EXAMPLE__
#define __TP_EXAMPLE__

#include "../ARTICULATED_RIGID_BODIES_3D_EXAMPLE.h"
namespace PhysBAM{

template<class T>
class TP_EXAMPLE:public ARTICULATED_RIGID_BODIES_3D_EXAMPLE<T,T>
{
public:
    RIGID_BODY<TV>* square1,*square2;

    TP_EXAMPLE(int parameter=0);
    ~TP_EXAMPLE();
    void Initialize_Rigid_Bodies();
    void Make_Lathe_Chain(VECTOR_3D<T>& start,VECTOR_3D<T> orient,int num_in_chain,bool start_joint,bool end_joint,int& num_joints,int& num_bodies);
    VECTOR_3D<T> Make_Lathe_Chain_2(VECTOR_3D<T> start_point,VECTOR_3D<T> direction,int number_of_links,bool start_joint,bool end_joint,int branch_joint,int branch_body,int branch_num,int& num_joints,int& num_bodies);
};
}
#endif
