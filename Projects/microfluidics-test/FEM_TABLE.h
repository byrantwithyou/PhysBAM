//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __FEM_TABLE__
#define __FEM_TABLE__
#include "COMMON.h"

namespace PhysBAM{

// Order[12]: u0.x, u1.x, u2.x, u3.x, u4.x, u5.x, u0.y, u1.y, u2.y, u3.y, u4.y, u5.y,
// 0-4, 1-5, 2-6 are opposing pairs; 0-2 are vertices, 3-5 are edges.
//
// Order[10]: u.x^2, v.x*u.x, u.y*u.x, u.x*v.y, v.x^2, u.y*v.x, v.y*v.x, u.y^2, v.y*u.y, v.y^2
// u = x1 - x0, v = x2 - x0

extern int force_table[6][6];

int Unique_Entries_Visc(LOCAL_V_CODE_ID i,int b);
int Unique_Entries_Pres(LOCAL_P_CODE_ID i,int b);
LOCAL_V_CODE_ID Number_Unique_Visc_Codes();
LOCAL_P_CODE_ID Number_Unique_Pres_Codes();

LOCAL_V_CODE_ID Main_Table_Visc(int dof_u,int dim_u,int dof_v,int dim_v);
LOCAL_V_CODE_ID Vertex_Table_Visc(int dim_u,int dim_v,int num_tri);
LOCAL_V_CODE_ID Edge_Table_Visc(int can_u,int can_v,int dof_u,int dim_u,int dof_v,int dim_v);

LOCAL_P_CODE_ID Main_Table_Pres(int dof_p,int dof_u,int dim_u);
LOCAL_P_CODE_ID Vertex_Table_Pres(int dim_u,int neg_x,int neg_z);
LOCAL_P_CODE_ID Edge_Table_Pres(int can_p,int can_u,int dof_p,int dof_u,int dim_u);

}
#endif
