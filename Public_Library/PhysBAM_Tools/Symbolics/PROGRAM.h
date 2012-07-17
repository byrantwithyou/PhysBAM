//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYREGHT.txt.
//#####################################################################
// Class PROGRAM
//#####################################################################
#ifndef __PROGRAM__
#define __PROGRAM__

#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <string>
namespace PhysBAM{

enum op_type {
    op_nop,op_copy,op_add,op_sub,op_mul,op_div,op_neg,op_inv,
    op_sqrt,op_exp,op_ln,op_pow,
    op_lt,op_le,op_gt,op_ge,op_eq,op_ne,op_not,op_or,op_and,
    op_br_z,op_br_nz,op_jmp,op_label,
    op_last
};

struct INSTRUCTION
{
    op_type type;
    int dest,src0,src1;
};

/*
  constants(0)=0
  constants(1)=1
 */
template<class T> struct PROGRAM_CONTEXT;
template<class T>
struct PROGRAM
{
    int extra_out;
    int num_tmp;
    int num_labels;
    ARRAY<T> constants;
    ARRAY<INSTRUCTION> code;
    ARRAY<std::string> var_in,var_out;
    HASHTABLE<std::string,int> dict;

    PROGRAM()
        :extra_out(0),num_tmp(0),num_labels(0)
    {
        constants.Append(0);
        constants.Append(1);
    }
    void Execute(ARRAY<T>& reg) const;
    void Execute(PROGRAM_CONTEXT<T>& context) const
    {Execute(context.reg);}
    void Execute_Op(ARRAY<T>& reg,int& ip) const;
    int Diff(int diff_expr,int diff_var);
    void Parse(const char* str);
    int Parse_Command(const char*& str);
    void Print() const;
    void Optimize();
    void Finalize();
};
template<class T>
struct PROGRAM_CONTEXT
{
    ARRAY<T> reg;
    ARRAY_VIEW<T> data_in,data_out;

    PROGRAM_CONTEXT(const PROGRAM<T>& prog)
        :reg(prog.num_tmp+prog.var_in.m+prog.var_out.m+prog.extra_out+prog.constants.m),
        data_in(reg.Array_View(prog.num_tmp,prog.var_in.m)),
        data_out(reg.Array_View(prog.num_tmp+prog.var_in.m,prog.var_out.m+prog.extra_out))
    {
        reg.Array_View(prog.num_tmp+prog.var_in.m+prog.var_out.m+prog.extra_out,prog.constants.m)=prog.constants;
    }
};

}
#endif
