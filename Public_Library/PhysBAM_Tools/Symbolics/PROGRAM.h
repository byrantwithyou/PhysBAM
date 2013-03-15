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
#include <map>
#include <string>
namespace PhysBAM{

enum op_type {
    op_nop,op_copy,op_add,op_sub,op_mul,op_div,op_neg,op_inv,
    op_sqrt,op_exp,op_ln,op_pow,
    op_lt,op_le,op_gt,op_ge,op_eq,op_ne,op_not,op_or,op_and,
    op_br_z,op_br_nz,op_jmp,op_label,op_phi,
    op_last
};

const int mem_shift=28;
const int mem_mask=0xF0000000;
const int mem_reg=0x00000000;
const int mem_const=0x10000000;
const int mem_in=0x20000000;
const int mem_out=0x30000000;

struct INSTRUCTION
{
    op_type type;
    int dest,src0,src1;
};

template<class T>
struct CODE_BLOCK
{
    CODE_BLOCK *prev[2],*next[2];
    ARRAY<INSTRUCTION> code;
    HASHTABLE<int> in,out;
    int id;

    CODE_BLOCK(){prev[0]=0;prev[1]=0;next[0]=0;next[1]=0;}
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
    std::map<T,int> constant_lookup;
    ARRAY<CODE_BLOCK<T>*> code_blocks;
    ARRAY<INSTRUCTION> flat_code;
    ARRAY<std::string> var_in,var_out;
    HASHTABLE<std::string,int> dict;
    bool finalized;
    VECTOR<ARRAY<VECTOR<int,2> >,4> defs;
    VECTOR<ARRAY<HASHTABLE<VECTOR<int,2> > >,4> uses;
    bool def_use_stale;

    struct VARIABLE_STATE
    {
        T value;
        enum state_type {unknown,constant,overdetermined} state;
    };

    PROGRAM()
        :extra_out(0),num_tmp(0),num_labels(0),finalized(false),def_use_stale(true)
    {
    }
    void Execute(ARRAY<T>& reg) const;
    void Execute(PROGRAM_CONTEXT<T>& context) const
    {Execute(context.reg);}
    T Evaluate_Op(int type,T in0,T in1) const;
    bool Deduce_Op(int type,T& out,T* in0,T* in1) const;
    void Execute_Op(ARRAY<T>& reg,int& ip) const;
    int Diff(int diff_expr,int diff_var);
    void Parse(const char* str,bool keep_all_vars=false);
    void Print() const;
    void Print(const ARRAY<INSTRUCTION>& code) const;
    void Optimize();
    int Add_Constant(T x);
    void Finalize();

    VECTOR<int,2> Get_Def(int var) const {return defs(var>>mem_shift)(var&~mem_mask);}
    void Set_Def(int var,const VECTOR<int,2>& def) {defs(var>>mem_shift)(var&~mem_mask)=def;}
    const HASHTABLE<VECTOR<int,2> >& Get_Uses(int var) const {return uses(var>>mem_shift)(var&~mem_mask);}
    HASHTABLE<VECTOR<int,2> >& Get_Uses(int var) {return uses(var>>mem_shift)(var&~mem_mask);}
    void Add_Use(int var,const VECTOR<int,2>& use) {uses(var>>mem_shift)(var&~mem_mask).Set(use);}
    const HASHTABLE<VECTOR<int,2> >& Get_Uses(const VECTOR<int,2>& I) const {return Get_Uses(Get_Instruction(I).dest);}
    INSTRUCTION& Get_Instruction(const VECTOR<int,2>& I) {return code_blocks(I.x)->code(I.y);}
    const INSTRUCTION& Get_Instruction(const VECTOR<int,2>& I) const {return code_blocks(I.x)->code(I.y);}

protected:
    int Parse_Command(const char*& str);
    void Make_SSA();
    void Make_SSA_Relabel(int& count,ARRAY<ARRAY<int> >& S,HASHTABLE<int>& V_used,const ARRAY<ARRAY<int> >& idom_tree_children,int bl);
    void Leave_SSA();
    void Compute_In_Out(CODE_BLOCK<T>* block);
    void Compute_In_Out();
    void Relabel_Registers(ARRAY<int>& var_map);
    void Update_Use_Def();
    void Copy_Propagation();
    void Simplify_Phis();
    void Sparse_Conditional_Constant_Propagation();
    void SCCP_Visit_Instruction(const VECTOR<int,2>& I,VECTOR<ARRAY<VARIABLE_STATE>,4>& variable_state,
        ARRAY<CODE_BLOCK<T>*>& block_worklist,ARRAY<VECTOR<int,2> >& op_worklist,ARRAY<bool>& block_exec);
    void Remove_Dead_Code();
    void Detach_Next(CODE_BLOCK<T>* B,int j);
    void Propagate_Copy(int old_var,int new_var);
    void Local_Common_Expresssion_Elimination(CODE_BLOCK<T>* B);
    void Local_Common_Expresssion_Elimination();
    void Reduce_In_Place(INSTRUCTION& o);
    void Compress_Registers();
    void Eliminate_Unused_Register(int var);
    void Eliminate_Unused_Constant(int var);
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
