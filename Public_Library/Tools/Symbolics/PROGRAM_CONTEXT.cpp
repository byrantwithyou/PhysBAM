//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYREGHT.txt.
//#####################################################################
#include <Tools/Symbolics/PROGRAM.h>
#include <Tools/Symbolics/PROGRAM_CONTEXT.h>
namespace PhysBAM{
//#####################################################################
// Function Initialize
//#####################################################################
template<class T> void PROGRAM_CONTEXT<T>::
Initialize(const PROGRAM<T>& prog)
{
    reg.Resize(prog.num_tmp+prog.var_in.m+prog.var_out.m+prog.extra_out+prog.constants.m);
    data_in.Set(reg.Array_View(prog.num_tmp,prog.var_in.m));
    data_out.Set(reg.Array_View(prog.num_tmp+prog.var_in.m,prog.var_out.m+prog.extra_out));
    reg.Array_View(prog.num_tmp+prog.var_in.m+prog.var_out.m+prog.extra_out,prog.constants.m)=prog.constants;
}
template struct PROGRAM_CONTEXT<float>;
template struct PROGRAM_CONTEXT<double>;
}
