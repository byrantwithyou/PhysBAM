//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __FEM_TABLE__
#define __FEM_TABLE__
#include <Core/Matrices/MATRIX.h>
#include "COMMON.h"

namespace PhysBAM{

// Order[12]: u0.x, u1.x, u2.x, u3.x, u4.x, u5.x, u0.y, u1.y, u2.y, u3.y, u4.y, u5.y,
// 0-4, 1-5, 2-6 are opposing pairs; 0-2 are vertices, 3-5 are edges.
//
// Order[10]: u.x^2, v.x*u.x, u.y*u.x, u.x*v.y, v.x^2, u.y*v.x, v.y*v.x, u.y^2, v.y*u.y, v.y^2
// u = x1 - x0, v = x2 - x0

int Unique_Entries_Visc(LOCAL_V_CODE_ID i,int b);
int Unique_Entries_Pres(LOCAL_P_CODE_ID i,int b);
LOCAL_V_CODE_ID Number_Unique_Visc_Codes();
LOCAL_P_CODE_ID Number_Unique_Pres_Codes();

MATRIX<LOCAL_V_CODE_ID,2> Main_Table_Visc(int dof_u,int dof_v);
MATRIX<LOCAL_V_CODE_ID,2> Vertex_Table_Visc(int num_tri);
MATRIX<LOCAL_V_CODE_ID,2> Edge_Table_Visc(int can_u,int can_v,int dof_u,int dof_v);

VECTOR<LOCAL_P_CODE_ID,2> Main_Table_Pres(int dof_p,int dof_u);
VECTOR<LOCAL_P_CODE_ID,2> Vertex_Table_Pres(int neg_x,int neg_z);
VECTOR<LOCAL_P_CODE_ID,2> Edge_Table_Pres(int can_p,int can_u,int dof_p,int dof_u);

template<class T,class TV>
VECTOR<TV,6> Times_force_NdotN(const VECTOR<TV,6>& f, T tri_area);

// Result should be scaled by edge length
template<class TV>
VECTOR<TV,3> Times_BC_NdotN(const VECTOR<TV,3>& bc);

template<class T>
VECTOR<T,3> Times_div_PdotN(const VECTOR<T,6>& div, T tri_area);

}
#endif
