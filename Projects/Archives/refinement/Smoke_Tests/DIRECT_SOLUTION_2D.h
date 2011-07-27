//#####################################################################
// Copyright 2009, Wen Zheng.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __DIRECT_SOLUTION_2D__
#define __DIRECT_SOLUTION_2D__

#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
namespace PhysBAM{

int Get_Boundary_Condition_Code(const ARRAY<bool,FACE_INDEX<2> >& psi_N)
{
    typedef VECTOR<int,2> TV_INT;
    int code=0;

    if(!psi_N(FACE_INDEX<2>(1,TV_INT(2,1)))) code+=1;
    if(!psi_N(FACE_INDEX<2>(1,TV_INT(2,2)))) code+=2;

    if(!psi_N(FACE_INDEX<2>(2,TV_INT(1,2)))) code+=4;
    if(!psi_N(FACE_INDEX<2>(2,TV_INT(2,2)))) code+=8;

    return code;
}

template<class T>
void Solve_For_Pressure_Analytically(ARRAY<T,VECTOR<int,2> >& p,const ARRAY<T,VECTOR<int,2> >& div,const int& code)
{
typedef VECTOR<int,2> TV_INT;

switch(code){
case 0:
p(TV_INT(1,1))=0;
p(TV_INT(1,2))=0;
p(TV_INT(2,1))=0;
p(TV_INT(2,2))=0;
break;

case 1:
p(TV_INT(1,1))=T(-1)*div(TV_INT(1,1));
p(TV_INT(1,2))=0;
p(TV_INT(2,1))=0;
p(TV_INT(2,2))=0;
break;

case 2:
p(TV_INT(1,1))=0;
p(TV_INT(1,2))=0;
p(TV_INT(2,1))=T(-1)*div(TV_INT(2,1));
p(TV_INT(2,2))=0;
break;

case 3:
p(TV_INT(1,1))=T(-1)*div(TV_INT(1,1));
p(TV_INT(1,2))=0;
p(TV_INT(2,1))=T(-1)*div(TV_INT(2,1));
p(TV_INT(2,2))=0;
break;

case 4:
p(TV_INT(1,1))=T(-1)*div(TV_INT(1,1));
p(TV_INT(1,2))=0;
p(TV_INT(2,1))=0;
p(TV_INT(2,2))=0;
break;

case 5:
p(TV_INT(1,1))=T(-1)*div(TV_INT(1,1))+T(-1)*div(TV_INT(1,2));
p(TV_INT(1,2))=T(-1)*div(TV_INT(1,1))+T(-2)*div(TV_INT(1,2));
p(TV_INT(2,1))=0;
p(TV_INT(2,2))=0;
break;

case 6:
p(TV_INT(1,1))=T(-2)*div(TV_INT(1,1))+T(-1)*div(TV_INT(2,1));
p(TV_INT(1,2))=0;
p(TV_INT(2,1))=T(-1)*div(TV_INT(1,1))+T(-1)*div(TV_INT(2,1));
p(TV_INT(2,2))=0;
break;

case 7:
p(TV_INT(1,1))=T(-2)*div(TV_INT(1,1))+T(-2)*div(TV_INT(1,2))+T(-1)*div(TV_INT(2,1));
p(TV_INT(1,2))=T(-2)*div(TV_INT(1,1))+T(-3)*div(TV_INT(1,2))+T(-1)*div(TV_INT(2,1));
p(TV_INT(2,1))=T(-1)*div(TV_INT(1,1))+T(-1)*div(TV_INT(1,2))+T(-1)*div(TV_INT(2,1));
p(TV_INT(2,2))=0;
break;

case 8:
p(TV_INT(1,1))=0;
p(TV_INT(1,2))=T(-1)*div(TV_INT(1,2));
p(TV_INT(2,1))=0;
p(TV_INT(2,2))=0;
break;

case 9:
p(TV_INT(1,1))=T(-2)*div(TV_INT(1,1))+T(-1)*div(TV_INT(1,2));
p(TV_INT(1,2))=T(-1)*div(TV_INT(1,1))+T(-1)*div(TV_INT(1,2));
p(TV_INT(2,1))=0;
p(TV_INT(2,2))=0;
break;

case 10:
p(TV_INT(1,1))=0;
p(TV_INT(1,2))=T(-1)*div(TV_INT(1,2));
p(TV_INT(2,1))=T(-1)*div(TV_INT(2,1));
p(TV_INT(2,2))=0;
break;

case 11:
p(TV_INT(1,1))=T(-2)*div(TV_INT(1,1))+T(-1)*div(TV_INT(1,2));
p(TV_INT(1,2))=T(-1)*div(TV_INT(1,1))+T(-1)*div(TV_INT(1,2));
p(TV_INT(2,1))=T(-1)*div(TV_INT(2,1));
p(TV_INT(2,2))=0;
break;

case 12:
p(TV_INT(1,1))=T(-1)*div(TV_INT(1,1));
p(TV_INT(1,2))=T(-1)*div(TV_INT(1,2));
p(TV_INT(2,1))=0;
p(TV_INT(2,2))=0;
break;

case 13:
p(TV_INT(1,1))=T(-2)*div(TV_INT(1,1))+T(-1)*div(TV_INT(1,2))+T(-2)*div(TV_INT(2,1));
p(TV_INT(1,2))=T(-1)*div(TV_INT(1,1))+T(-1)*div(TV_INT(1,2))+T(-1)*div(TV_INT(2,1));
p(TV_INT(2,1))=T(-2)*div(TV_INT(1,1))+T(-1)*div(TV_INT(1,2))+T(-3)*div(TV_INT(2,1));
p(TV_INT(2,2))=0;
break;

case 14:
p(TV_INT(1,1))=T(-2)*div(TV_INT(1,1))+T(-1)*div(TV_INT(2,1));
p(TV_INT(1,2))=T(-1)*div(TV_INT(1,2));
p(TV_INT(2,1))=T(-1)*div(TV_INT(1,1))+T(-1)*div(TV_INT(2,1));
p(TV_INT(2,2))=0;
break;

case 15:
p(TV_INT(1,1))=T(-1)*div(TV_INT(1,1))+T(-1./2.)*div(TV_INT(1,2))+T(-1./2.)*div(TV_INT(2,1));
p(TV_INT(1,2))=T(-1./2.)*div(TV_INT(1,1))+T(-3./4.)*div(TV_INT(1,2))+T(-1./4.)*div(TV_INT(2,1));
p(TV_INT(2,1))=T(-1./2.)*div(TV_INT(1,1))+T(-1./4.)*div(TV_INT(1,2))+T(-3./4.)*div(TV_INT(2,1));
p(TV_INT(2,2))=0;
break;
}
}
}
#endif
