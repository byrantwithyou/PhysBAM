//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYREGHT.txt.
//#####################################################################
// Class PROGRAM
//#####################################################################
#ifndef __PROGRAM__
#define __PROGRAM__

#include <Tools/Data_Structures/HASHTABLE.h>
#include <Tools/Symbolics/CODE_BLOCK.h>
#include <Tools/Symbolics/CODE_BLOCK_NODE.h>
#include <Tools/Symbolics/INSTRUCTION.h>
#include <Tools/Symbolics/PROGRAM_CONTEXT.h>
#include <Tools/Symbolics/PROGRAM_DEFINITIONS.h>
#include <Tools/Utilities/NONCOPYABLE.h>
#include <Tools/Vectors/VECTOR.h>
#include <map>
#include <string>
namespace PhysBAM{

template<class T>
struct DISTRIBUTE_HELPER
{
    T coeff;
    ARRAY<int> product;
};

/*
  constants(0)=0
  constants(1)=1
 */
template<class T>
struct PROGRAM
{
    int extra_out;
    int num_tmp;
    int num_labels;
    int last_block;
    ARRAY<T> constants;
    std::map<T,int> constant_lookup;
    ARRAY<CODE_BLOCK*> code_blocks;
    ARRAY<INSTRUCTION> flat_code;
    ARRAY<std::string> var_in,var_out;
    HASHTABLE<std::string,int> dict;
    bool finalized;
    VECTOR<ARRAY<CODE_BLOCK_NODE*>,4> defs;
    VECTOR<ARRAY<HASHTABLE<CODE_BLOCK_NODE*> >,4> uses;
    bool debug;

    struct VARIABLE_STATE
    {
        T value;
        enum state_type {unknown,constant,overdetermined} state;
    };

    PROGRAM()
        :extra_out(0),num_tmp(0),num_labels(0),last_block(0),finalized(false),debug(false)
    {
    }
    void Execute(ARRAY<T>& reg) const;
    void Execute(PROGRAM_CONTEXT<T>& context) const
    {Execute(context.reg);}
    T Evaluate_Op(int type,T in0,T in1) const;
    bool Deduce_Op(int type,T& out,T* in0,T* in1) const;
    void Execute_Op(ARRAY<T>& reg,int& ip) const;
    int Diff(const ARRAY<int>& diff_expr,int diff_var);
    void Parse(const char* str,bool keep_all_vars=false);
    void Print() const;
    void Print(const INSTRUCTION& I) const;
    void Print(const ARRAY<INSTRUCTION>& code) const;
    void Print(const CODE_BLOCK* B) const;
    void Optimize();
    int Add_Constant(T x);
    void Finalize();

    CODE_BLOCK_NODE* Get_Def(int var) const {return defs(var>>mem_shift)(var&~mem_mask);}
    void Set_Def(int var,CODE_BLOCK_NODE* def) {defs(var>>mem_shift)(var&~mem_mask)=def;}
    const HASHTABLE<CODE_BLOCK_NODE*>& Get_Uses(int var) const {return uses(var>>mem_shift)(var&~mem_mask);}
    HASHTABLE<CODE_BLOCK_NODE*>& Get_Uses(int var) {return uses(var>>mem_shift)(var&~mem_mask);}
    const HASHTABLE<CODE_BLOCK_NODE*>& Get_Uses(const CODE_BLOCK_NODE* N) const {return Get_Uses(N->inst.dest);}
    void Remove_Use_Def(CODE_BLOCK_NODE* N);
    void Add_Use_Def(CODE_BLOCK_NODE* N);
    void Delete_Instruction(CODE_BLOCK_NODE* N);

protected:
    void Parse_Command(const char* str);
    int Process_Node(PROGRAM_PARSE_NODE* node);
    void Print_Node(PROGRAM_PARSE_NODE* node);

    void Make_SSA();
    void Make_SSA_Relabel(int& count,ARRAY<ARRAY<int> >& S,HASHTABLE<int>& V_used,const ARRAY<ARRAY<int> >& idom_tree_children,int bl);
    void Leave_SSA();
    void Relabel_Registers(ARRAY<int>& var_map);
    void Update_Use_Def();
    void Copy_Propagation();
    void Simplify_Phis();
    void Sparse_Conditional_Constant_Propagation();
    void SCCP_Visit_Instruction(CODE_BLOCK_NODE* N,VECTOR<ARRAY<VARIABLE_STATE>,4>& variable_state,
        ARRAY<CODE_BLOCK*>& block_worklist,ARRAY<CODE_BLOCK_NODE*>& op_worklist,ARRAY<bool>& block_exec);
    void Remove_Dead_Code();
    void Detach_Next(CODE_BLOCK* B,int j);
    void Propagate_Copy(int old_var,int new_var);
    void Local_Common_Expresssion_Elimination(CODE_BLOCK* B);
    void Local_Common_Expresssion_Elimination();
    void Reduce_In_Place(CODE_BLOCK_NODE* N);
    void Compress_Registers();
    void Eliminate_Unused_Register(int var);
    void Eliminate_Unused_Constant(int var);
    CODE_BLOCK_NODE* Add_Raw_Instruction_To_Block_After(CODE_BLOCK* B,CODE_BLOCK_NODE* N,op_type type,int dest,int src0,int src1);
    CODE_BLOCK_NODE* Add_Raw_Instruction_To_Block_Before(CODE_BLOCK* B,CODE_BLOCK_NODE* N,op_type type,int dest,int src0,int src1);
    int Append_Instruction(op_type type,int dest,int src0,int src1);
    void Propagate_Copy_And_Remove_Instruction(CODE_BLOCK_NODE* N);
    void Simplify_With_Distributive_Law_Helper_Prod(int var,ARRAY<CODE_BLOCK_NODE*>& extra_nodes,DISTRIBUTE_HELPER<T>& prod);
    void Simplify_With_Distributive_Law_Helper(int var,T coeff,ARRAY<CODE_BLOCK_NODE*>& extra_nodes,ARRAY<DISTRIBUTE_HELPER<T> >& summands);
    int Simplify_With_Distributive_Law_Emit_Prod(DISTRIBUTE_HELPER<T>& prod,CODE_BLOCK_NODE* before,bool& need_neg);
    int Simplify_With_Distributive_Law_Emit_Simple(ARRAY<DISTRIBUTE_HELPER<T> >& summands,CODE_BLOCK_NODE* before,bool& need_neg);
    int Simplify_With_Distributive_Law_Emit_Consts(ARRAY<DISTRIBUTE_HELPER<T> >& summands,CODE_BLOCK_NODE* before,bool& need_neg);
    int Simplify_With_Distributive_Law_Emit(ARRAY<DISTRIBUTE_HELPER<T> >& summands,CODE_BLOCK_NODE* before,bool& need_neg);
    CODE_BLOCK_NODE* Simplify_With_Distributive_Law(CODE_BLOCK_NODE* N);
    void Simplify_With_Distributive_Law();
    void Change_Instruction(CODE_BLOCK_NODE* N,op_type type,int dest,int src0,int src1);
};
}
#endif
