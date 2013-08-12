#include <Tools/Arrays/INDIRECT_ARRAY.h>
#include <Tools/Data_Structures/HASHTABLE.h>
#include <Tools/Data_Structures/HASHTABLE_ITERATOR.h>
#include <Tools/Data_Structures/UNION_FIND.h>
#include <Tools/Log/LOG.h>
#include <Tools/Symbolics/CODE_BLOCK.h>
#include <Tools/Symbolics/DOMINATORS.h>
#include <Tools/Symbolics/INSTRUCTION.h>
#include <Tools/Symbolics/PROGRAM.h>
#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <set>
namespace PhysBAM{

enum op_flags {
    flag_none=0x00,
    flag_reg_dest=0x01,flag_reg_src0=0x02,flag_reg_src1=0x04,flag_can_deduce=0x08,
    flag_commute=0x10};

int op_flags_table[op_last];
int Init_Instructions()
{
    op_flags_table[op_nop]=flag_none;
    op_flags_table[op_copy]=flag_reg_dest|flag_reg_src0;
    op_flags_table[op_add]=flag_reg_dest|flag_reg_src0|flag_reg_src1|flag_commute;
    op_flags_table[op_sub]=flag_reg_dest|flag_reg_src0|flag_reg_src1;
    op_flags_table[op_mul]=flag_reg_dest|flag_reg_src0|flag_reg_src1|flag_can_deduce|flag_commute;
    op_flags_table[op_div]=flag_reg_dest|flag_reg_src0|flag_reg_src1|flag_can_deduce;
    op_flags_table[op_neg]=flag_reg_dest|flag_reg_src0;
    op_flags_table[op_inv]=flag_reg_dest|flag_reg_src0;
    op_flags_table[op_sqrt]=flag_reg_dest|flag_reg_src0;
    op_flags_table[op_exp]=flag_reg_dest|flag_reg_src0;
    op_flags_table[op_ln]=flag_reg_dest|flag_reg_src0;
    op_flags_table[op_pow]=flag_reg_dest|flag_reg_src0|flag_reg_src1|flag_can_deduce;
    op_flags_table[op_lt]=flag_reg_dest|flag_reg_src0|flag_reg_src1;
    op_flags_table[op_le]=flag_reg_dest|flag_reg_src0|flag_reg_src1;
    op_flags_table[op_gt]=flag_reg_dest|flag_reg_src0|flag_reg_src1;
    op_flags_table[op_ge]=flag_reg_dest|flag_reg_src0|flag_reg_src1;
    op_flags_table[op_eq]=flag_reg_dest|flag_reg_src0|flag_reg_src1|flag_commute;
    op_flags_table[op_ne]=flag_reg_dest|flag_reg_src0|flag_reg_src1|flag_commute;
    op_flags_table[op_not]=flag_reg_dest|flag_reg_src0;
    op_flags_table[op_or]=flag_reg_dest|flag_reg_src0|flag_reg_src1|flag_can_deduce|flag_commute;
    op_flags_table[op_and]=flag_reg_dest|flag_reg_src0|flag_reg_src1|flag_can_deduce|flag_commute;
    op_flags_table[op_br_z]=flag_reg_src0;
    op_flags_table[op_br_nz]=flag_reg_src0;
    op_flags_table[op_jmp]=flag_none;
    op_flags_table[op_label]=flag_none;
    op_flags_table[op_phi]=flag_reg_dest|flag_reg_src0|flag_reg_src1;
    return 0;
}
int call_init_instructions=Init_Instructions();
//#####################################################################
// Function Execute
//#####################################################################
template<class T> T PROGRAM<T>::
Evaluate_Op(int type,T in0,T in1) const
{
    switch(type){
        case op_copy: return in0;
        case op_add: return in0+in1;
        case op_sub: return in0-in1;
        case op_mul: return in0*in1;
        case op_div: return in0/in1;
        case op_neg: return -in0;
        case op_inv: return 1/in0;
        case op_sqrt: return sqrt(in0);
        case op_exp: return exp(in0);
        case op_ln: return log(in0);
        case op_pow: return pow(in0,in1);
        case op_lt: return in0<in1;
        case op_le: return in0<=in1;
        case op_gt: return in0>in1;
        case op_ge: return in0>=in1;
        case op_eq: return in0==in1;
        case op_ne: return in0!=in1;
        case op_not: return !in0;
        case op_or: return in0||in1;
        case op_and: return in0&&in1;
        default: PHYSBAM_FATAL_ERROR("Cannot evaluate this instruction");}
}
//#####################################################################
// Function Execute
//#####################################################################
template<class T> void PROGRAM<T>::
Execute_Op(ARRAY<T>& reg,int& ip) const
{
    const INSTRUCTION& o=flat_code(ip);
    if(op_flags_table[o.type]&flag_reg_dest){
        T in0=0,in1=0;
        if(op_flags_table[o.type]&flag_reg_src0) in0=reg(o.src0);
        if(op_flags_table[o.type]&flag_reg_src1) in1=reg(o.src1);
        reg(o.dest)=Evaluate_Op(o.type,in0,in1);
        return;}

    switch(o.type){
        case op_nop: case op_label: break;
        case op_br_z: if(!reg(o.src0)) ip=o.dest-1;break;
        case op_br_nz: if(reg(o.src0)) ip=o.dest-1;break;
        case op_jmp: ip=o.dest-1;break;
        default: PHYSBAM_FATAL_ERROR("Missing instruction");}
}
//#####################################################################
// Function Execute
//#####################################################################
template<class T> void PROGRAM<T>::
Execute(ARRAY<T>& reg) const
{
    PHYSBAM_ASSERT(finalized);
    for(int ip=0;ip<flat_code.m;ip++)
        Execute_Op(reg,ip);
}
//#####################################################################
// Function Evaluate_Op
//#####################################################################
template<class T> bool PROGRAM<T>::
Deduce_Op(int type,T& out,T* in0,T* in1) const
{
    switch(type){
        case op_mul: if((in0 && !*in0) || (in1 && !*in1)){out=0;return true;} return false;
        case op_div: if((in0 && !*in0)){out=0;return true;} return false;
        case op_pow: if(in0 && (*in0==0 || *in0==1)){out=*in0;return true;}
            if(in1 && *in1==0){out=1;return true;} return false;
        case op_or: if((in0 && *in0) || (in1 && *in1)){out=1;return true;} return false;
        case op_and: if((in0 && !*in0) || (in1 && !*in1)){out=0;return true;}} return false;
    return false;
}
//#####################################################################
// Function Diff
//#####################################################################
template<class T> CODE_BLOCK_NODE* PROGRAM<T>::
Add_Raw_Instruction_To_Block_Before(CODE_BLOCK* B,CODE_BLOCK_NODE* N,op_type type,int dest,int src0,int src1)
{
    INSTRUCTION in={type,dest,src0,src1};
    CODE_BLOCK_NODE* M=B->Insert_Before(in,N);
    if(op_flags_table[type]&flag_reg_dest && (dest&~mem_mask)>=uses(dest>>mem_shift).m)
        uses(dest>>mem_shift).Resize((dest&~mem_mask)+1);
    Add_Use_Def(M);
    return M;
}
//#####################################################################
// Function Diff
//#####################################################################
template<class T> CODE_BLOCK_NODE* PROGRAM<T>::
Add_Raw_Instruction_To_Block_After(CODE_BLOCK* B,CODE_BLOCK_NODE* N,op_type type,int dest,int src0,int src1)
{
    INSTRUCTION in={type,dest,src0,src1};
    CODE_BLOCK_NODE* M=B->Insert_After(in,N);
    if(op_flags_table[type]&flag_reg_dest && (dest&~mem_mask)>=uses(dest>>mem_shift).m)
        uses(dest>>mem_shift).Resize((dest&~mem_mask)+1);
    Add_Use_Def(M);
    return M;
}
//#####################################################################
// Function Diff
//#####################################################################
template<class T> int PROGRAM<T>::
Diff(int diff_expr,int diff_var)
{
    Optimize();

    ARRAY<int> diff_tmp(num_tmp);
    diff_tmp.Fill(-1);
    uses(3).Append(HASHTABLE<CODE_BLOCK_NODE*>());

    for(int bl=0;bl<code_blocks.m;bl++){
        CODE_BLOCK* B=code_blocks(bl);
        if(!B) continue;
        for(CODE_BLOCK_NODE* N=B->head;N;N=N->next){
            const INSTRUCTION& o=N->inst;
            if(!(op_flags_table[o.type]&flag_reg_dest)) continue;

            int d=-1,s0=-1,s1=-1;
            if((o.dest&mem_mask)==mem_out){
                if(o.dest!=(mem_out|diff_expr)) continue;
                d=(var_out.m+extra_out)|mem_out;}
            else if(diff_tmp(o.dest)>=0) d=diff_tmp(o.dest);
            else diff_tmp(o.dest)=d=num_tmp++;

            if(op_flags_table[o.type]&flag_reg_src0){
                if((o.src0&mem_mask)==mem_const) s0=Add_Constant(0);
                else if((o.src0&mem_mask)==mem_in) s0=Add_Constant(o.src0==(mem_in|diff_var));
                else s0=diff_tmp(o.src0);}

            if(op_flags_table[o.type]&flag_reg_src1){
                if((o.src1&mem_mask)==mem_const) s1=Add_Constant(0);
                else if((o.src1&mem_mask)==mem_in) s1=Add_Constant(o.src1==(mem_in|diff_var));
                else s1=diff_tmp(o.src1);}

            switch(o.type){
                case op_copy:case op_neg:N=Add_Raw_Instruction_To_Block_After(B,N,o.type,d,s0,-1);break;
                case op_add:case op_sub:case op_phi:N=Add_Raw_Instruction_To_Block_After(B,N,o.type,d,s0,s1);break;
                case op_mul:
                    N=Add_Raw_Instruction_To_Block_After(B,N,op_mul,num_tmp,s0,o.src1);
                    N=Add_Raw_Instruction_To_Block_After(B,N,op_mul,num_tmp+1,o.src0,s1);
                    N=Add_Raw_Instruction_To_Block_After(B,N,op_add,d,num_tmp,num_tmp+1);
                    num_tmp+=2;
                    break;
                case op_div:
                    N=Add_Raw_Instruction_To_Block_After(B,N,op_div,num_tmp,s1,o.src1);
                    N=Add_Raw_Instruction_To_Block_After(B,N,op_mul,num_tmp+1,o.dest,num_tmp);
                    N=Add_Raw_Instruction_To_Block_After(B,N,op_div,num_tmp+2,s0,o.src1);
                    N=Add_Raw_Instruction_To_Block_After(B,N,op_sub,d,num_tmp+2,num_tmp+1);
                    num_tmp+=3;
                    break;
                case op_inv:
                    N=Add_Raw_Instruction_To_Block_After(B,N,op_mul,num_tmp,o.dest,o.dest);
                    N=Add_Raw_Instruction_To_Block_After(B,N,op_mul,num_tmp+1,s0,num_tmp);
                    N=Add_Raw_Instruction_To_Block_After(B,N,op_neg,d,num_tmp+1,-1);
                    num_tmp+=2;
                    break;
                case op_sqrt:
                    N=Add_Raw_Instruction_To_Block_After(B,N,op_div,num_tmp,s0,o.dest);
                    N=Add_Raw_Instruction_To_Block_After(B,N,op_mul,d,Add_Constant((T).5),num_tmp++);
                    break;
                case op_exp:N=Add_Raw_Instruction_To_Block_After(B,N,op_mul,d,s0,o.dest);break;
                case op_ln:N=Add_Raw_Instruction_To_Block_After(B,N,op_div,d,s0,o.src0);break;
                case op_pow:
                    N=Add_Raw_Instruction_To_Block_After(B,N,op_ln,num_tmp,o.src0,-1);
                    N=Add_Raw_Instruction_To_Block_After(B,N,op_mul,num_tmp+1,s1,num_tmp);
                    N=Add_Raw_Instruction_To_Block_After(B,N,op_mul,num_tmp+2,o.dest,num_tmp+1);
                    N=Add_Raw_Instruction_To_Block_After(B,N,op_sub,num_tmp+3,o.src1,Add_Constant(1));
                    N=Add_Raw_Instruction_To_Block_After(B,N,op_pow,num_tmp+4,o.src0,num_tmp+3);
                    N=Add_Raw_Instruction_To_Block_After(B,N,op_mul,num_tmp+5,num_tmp+4,o.src1);
                    N=Add_Raw_Instruction_To_Block_After(B,N,op_mul,num_tmp+6,num_tmp+5,s0);
                    N=Add_Raw_Instruction_To_Block_After(B,N,op_add,d,num_tmp+2,num_tmp+6);
                    num_tmp+=7;
                    break;
                case op_lt:case op_le:case op_gt:case op_ge:case op_eq:case op_ne:case op_not:case op_or:case op_and:
                    N=Add_Raw_Instruction_To_Block_After(B,N,op_copy,d,Add_Constant(0),-1);
                    break;
                default: PHYSBAM_FATAL_ERROR("Missing diff instruction");
            }
        }
    }

    if(debug) LOG::cout<<__FUNCTION__<<std::endl;
    Print();

    return var_out.m+extra_out++;
}

enum op_parse_info_flags {flag_dest_is_src0=1,flag_num_ok_after=2};

const char * op_chars="+-!*/^<>=!|&";

struct OPERATOR_DEFINITIONS
{
    HASHTABLE<std::string,OP_PARSE_INFO> op_lookup[2];
    HASHTABLE<std::string,op_type> func_lookup;

    bool right_to_left[prec_hi+1];
    
    OPERATOR_DEFINITIONS()
    {
        op_lookup[1].Set("+",OP_PARSE_INFO(op_copy,prec_neg,flag_num_ok_after));
        op_lookup[1].Set("-",OP_PARSE_INFO(op_neg,prec_neg,flag_num_ok_after));
        op_lookup[1].Set("!",OP_PARSE_INFO(op_not,prec_neg,flag_num_ok_after));

        op_lookup[0].Set("+",OP_PARSE_INFO(op_add,prec_add,flag_num_ok_after));
        op_lookup[0].Set("-",OP_PARSE_INFO(op_sub,prec_add,flag_num_ok_after));
        op_lookup[0].Set("*",OP_PARSE_INFO(op_mul,prec_mul,flag_num_ok_after));
        op_lookup[0].Set("/",OP_PARSE_INFO(op_div,prec_mul,flag_num_ok_after));
        op_lookup[0].Set("^",OP_PARSE_INFO(op_pow,prec_neg,flag_num_ok_after));
        op_lookup[0].Set("<",OP_PARSE_INFO(op_lt,prec_comp,flag_num_ok_after));
        op_lookup[0].Set("<=",OP_PARSE_INFO(op_le,prec_comp,flag_num_ok_after));
        op_lookup[0].Set(">",OP_PARSE_INFO(op_gt,prec_comp,flag_num_ok_after));
        op_lookup[0].Set(">=",OP_PARSE_INFO(op_ge,prec_comp,flag_num_ok_after));
        op_lookup[0].Set("==",OP_PARSE_INFO(op_eq,prec_eq,flag_num_ok_after));
        op_lookup[0].Set("!=",OP_PARSE_INFO(op_ne,prec_eq,flag_num_ok_after));
        op_lookup[0].Set("||",OP_PARSE_INFO(op_or,prec_or,flag_num_ok_after));
        op_lookup[0].Set("&&",OP_PARSE_INFO(op_and,prec_and,flag_num_ok_after));
        op_lookup[0].Set("=",OP_PARSE_INFO(op_copy,prec_assign,flag_num_ok_after));
        op_lookup[0].Set("+=",OP_PARSE_INFO(op_add,prec_assign,flag_dest_is_src0|flag_num_ok_after));
        op_lookup[0].Set("-=",OP_PARSE_INFO(op_sub,prec_assign,flag_dest_is_src0|flag_num_ok_after));
        op_lookup[0].Set("*=",OP_PARSE_INFO(op_mul,prec_assign,flag_dest_is_src0|flag_num_ok_after));
        op_lookup[0].Set("/=",OP_PARSE_INFO(op_div,prec_assign,flag_dest_is_src0|flag_num_ok_after));
        op_lookup[0].Set("^=",OP_PARSE_INFO(op_pow,prec_assign,flag_dest_is_src0|flag_num_ok_after));

        func_lookup.Set("inv",op_inv);
        func_lookup.Set("sqrt",op_sqrt);
        func_lookup.Set("exp",op_exp);
        func_lookup.Set("ln",op_ln);

        right_to_left[prec_lo]=false;
        right_to_left[prec_assign]=true;
        right_to_left[prec_cond]=true;
        right_to_left[prec_or]=false;
        right_to_left[prec_and]=false;
        right_to_left[prec_eq]=false;
        right_to_left[prec_comp]=false;
        right_to_left[prec_add]=false;
        right_to_left[prec_mul]=false;
        right_to_left[prec_neg]=true;
        right_to_left[prec_func]=false;
        right_to_left[prec_hi]=false;
    }
} op_tokens;
//#####################################################################
// Function Parse
//#####################################################################
template<class T> void PROGRAM<T>::
Parse(const char* str,bool keep_all_vars)
{
    code_blocks.Delete_Pointers_And_Clean_Memory();
    flat_code.Remove_All();
    finalized=false;

    CODE_BLOCK* A=new CODE_BLOCK;
    A->id=code_blocks.Append(A);

    for(int i=0;i<var_in.m;i++){
        int x=num_tmp++;
        INSTRUCTION in={op_copy,x,i|mem_in,-1};
        code_blocks(0)->Append(in);
        dict.Set(var_in(i),x);}
    for(int i=0;i<var_out.m;i++)
        if(!dict.Contains(var_out(i)))
            dict.Set(var_out(i),num_tmp++);

    if(debug) LOG::cout<<dict<<"  "<<num_tmp<<std::endl;

    Parse_Command(str);

    if(keep_all_vars){
        HASHTABLE<std::string> known;
        known.Set_All(var_out);
        ARRAY<std::string> extra_vars;
        for(HASHTABLE<std::string,int>::ITERATOR it(dict);it.Valid();it.Next())
            if(!known.Contains(it.Key()))
                extra_vars.Append(it.Key());
        extra_vars.Sort();
        var_out.Append_Elements(extra_vars);}

    if(debug) LOG::cout<<dict<<"  "<<num_tmp<<std::endl;

    if(debug) LOG::cout<<var_out<<std::endl;
    ARRAY<INSTRUCTION> zero_inst;
    for(int i=0;i<var_out.m;i++){
        INSTRUCTION in1={op_copy,i|mem_out,dict.Get(var_out(i)),-1};
        code_blocks(last_block)->Append(in1);
        INSTRUCTION in2={op_copy,dict.Get(var_out(i)),Add_Constant(0),-1};
        code_blocks(0)->Prepend(in2);}

    Print();

    Make_SSA();
}
//#####################################################################
// Function Pop_Op_Stack
//#####################################################################
template<class T> void PROGRAM<T>::
Pop_Op_Stack(ARRAY<CODE_BLOCK*>& block_stack,ARRAY<int>& data_stack,ARRAY<OP_PARSE_INFO>& op_stack)
{
    PHYSBAM_ASSERT(op_stack.m>0);
    OP_PARSE_INFO op=op_stack.Pop();
    if(op.prec==prec_cond){
        PHYSBAM_ASSERT(block_stack.m>0);
        CODE_BLOCK* C=block_stack.Pop();
        PHYSBAM_ASSERT(data_stack.m>1);
        int tos=data_stack.Pop();
        INSTRUCTION in_C={op_copy,data_stack.Last(),tos,-1};
        C->Append(in_C);
        return;}
    if(op.prec==prec_nest){LOG::cout<<"Expected '"<<(char)op.op<<"' but found ';'."<<std::endl;PHYSBAM_FATAL_ERROR("");}

    INSTRUCTION in={op.op,-1,-1,-1};
    if(op_flags_table[op.op]&flag_reg_src1){PHYSBAM_ASSERT(data_stack.m>0);in.src1=data_stack.Pop();}
    if(op_flags_table[op.op]&flag_reg_src0){PHYSBAM_ASSERT(data_stack.m>0);in.src0=data_stack.Pop();}
    if(op.op==op_copy && op.prec==prec_assign){PHYSBAM_ASSERT(data_stack.m>0);in.dest=data_stack.Pop();}
    else if(op.flags&flag_dest_is_src0) in.dest=in.src0;
    else if(op_flags_table[op.op]&flag_reg_dest) in.dest=num_tmp++;
    if(op_flags_table[op.op]&flag_reg_dest) data_stack.Append(in.dest);
    block_stack.Last()->Append(in);
}
//#####################################################################
// Function Unwind_Stack
//#####################################################################
template<class T> void PROGRAM<T>::
Unwind_Match(ARRAY<CODE_BLOCK*>& block_stack,ARRAY<int>& data_stack,ARRAY<OP_PARSE_INFO>& op_stack,char match)
{
    while(op_stack.m && op_stack.Last().prec>prec_nest) Pop_Op_Stack(block_stack,data_stack,op_stack);
    if(!op_stack.m){LOG::cout<<"Expected ';' but found '"<<match<<"'."<<std::endl;PHYSBAM_FATAL_ERROR("");}
    else if(op_stack.Last().prec!=prec_nest){LOG::cout<<"Invalid token on op stack."<<std::endl;PHYSBAM_FATAL_ERROR("");}
    else if(op_stack.Last().op!=match){LOG::cout<<"Expected '"<<(char)op_stack.Last().op<<"' but found '"<<match<<"'."<<std::endl;PHYSBAM_FATAL_ERROR("");}
    else op_stack.Pop();
}
//#####################################################################
// Function Unwind_Stack
//#####################################################################
template<class T> void PROGRAM<T>::
Unwind_To_Prec(ARRAY<CODE_BLOCK*>& block_stack,ARRAY<int>& data_stack,ARRAY<OP_PARSE_INFO>& op_stack,int prec)
{
    while(op_stack.m && op_stack.Last().prec>prec)
        Pop_Op_Stack(block_stack,data_stack,op_stack);
    if(!op_tokens.right_to_left[prec]) 
        while(op_stack.m && op_stack.Last().prec==prec)
            Pop_Op_Stack(block_stack,data_stack,op_stack);
}
//#####################################################################
// Function Unwind_Stack
//#####################################################################
template<class T> void PROGRAM<T>::
Unwind_Op_Stack(ARRAY<CODE_BLOCK*>& block_stack,ARRAY<int>& data_stack,ARRAY<OP_PARSE_INFO>& op_stack)
{
    while(op_stack.m)
        Pop_Op_Stack(block_stack,data_stack,op_stack);
}
//#####################################################################
// Function Unwind_Stack
//#####################################################################
template<class T> bool PROGRAM<T>::
Try_Operator(ARRAY<CODE_BLOCK*>& block_stack,ARRAY<int>& data_stack,ARRAY<OP_PARSE_INFO>& op_stack,bool& num_ok,const std::string& str)
{
    OP_PARSE_INFO opi;
    if(op_tokens.op_lookup[num_ok].Get(str,opi)){
        Unwind_To_Prec(block_stack,data_stack,op_stack,opi.prec);
        num_ok=opi.flags&flag_num_ok_after;
        op_stack.Append(opi);
        return true;}

    if(op_tokens.op_lookup[1-num_ok].Contains(str)){
        LOG::cout<<"Unexpected operator '"<<str<<"'."<<std::endl;
        PHYSBAM_FATAL_ERROR();}
    return false;
}
//#####################################################################
// Function Parse_Command
//#####################################################################
template<class T> void PROGRAM<T>::
Parse_Command(const char* str)
{
    bool num_ok=true;
    ARRAY<int> data_stack;
    ARRAY<OP_PARSE_INFO> op_stack;
    ARRAY<CODE_BLOCK*> block_stack;

    block_stack.Append(code_blocks.Last());

    while(1)
    {
        if(*str==';' || !*str){
            Unwind_Op_Stack(block_stack,data_stack,op_stack);

            if(data_stack.m>1){
                LOG::cout<<"Detected extra tokens on data stack while parsing."<<std::endl;
                PHYSBAM_FATAL_ERROR();}

            if(block_stack.m!=1){
                LOG::cout<<"Detected extra tokens on block stack while parsing."<<std::endl;
                PHYSBAM_FATAL_ERROR();}

            if(num_ok){
                if(data_stack.m){
                    LOG::cout<<"Unexpected ';'."<<std::endl;
                    PHYSBAM_FATAL_ERROR();}
                if(!*str) break;
                str++;
                num_ok=true;
                continue;}

            if(!data_stack.m){
                LOG::cout<<"Empty data stack at end of statement."<<std::endl;
                PHYSBAM_FATAL_ERROR();}

            data_stack.Pop();
            if(!*str) break;
            str++;
            num_ok=true;
            continue;}

        if(isspace(*str)){str++;continue;}

        if(*str=='('){
            str++;
            OP_PARSE_INFO opi(')');
            op_stack.Append(opi);
            num_ok=true;
            continue;}

        if(*str==')'){
            if(num_ok){
                LOG::cout<<"Unexpected ')'."<<std::endl;
                PHYSBAM_FATAL_ERROR();}
            str++;
            Unwind_Match(block_stack,data_stack,op_stack,')');
            num_ok=false;
            continue;}

        if(*str=='?'){
            str++;
            if(num_ok){
                LOG::cout<<"Unexpected '?'."<<std::endl;
                PHYSBAM_FATAL_ERROR();}
            Unwind_To_Prec(block_stack,data_stack,op_stack,prec_cond);
            num_ok=true;

            CODE_BLOCK* A=block_stack.Pop();
            CODE_BLOCK* B=new CODE_BLOCK;
            CODE_BLOCK* C=new CODE_BLOCK;
            CODE_BLOCK* D=new CODE_BLOCK;
            for(int s=0;s<2;s++){
                CODE_BLOCK* P=A->next[s];
                D->next[s]=P;
                if(P) P->prev[P->prev[1]==A]=D;}

            A->next[0]=B;
            A->next[1]=C;
            D->prev[0]=B;
            D->prev[1]=C;
            B->prev[0]=A;
            B->next[0]=D;
            C->prev[0]=A;
            C->next[0]=D;
            B->id=code_blocks.Append(B);
            C->id=code_blocks.Append(C);
            D->id=code_blocks.Append(D);
            block_stack.Append(D);
            block_stack.Append(C);
            block_stack.Append(B);

            INSTRUCTION in_A={op_br_z,C->id,data_stack.Pop(),B->id};
            A->Append(in_A);
            data_stack.Append(num_tmp++);

            op_stack.Append(OP_PARSE_INFO(':'));
            continue;}

        if(*str==':'){
            str++;
            if(num_ok){
                LOG::cout<<"Unexpected ':'."<<std::endl;
                PHYSBAM_FATAL_ERROR();}
            Unwind_Match(block_stack,data_stack,op_stack,':');
            num_ok=true;

            CODE_BLOCK* B=block_stack.Pop();
            int tos=data_stack.Pop();
            INSTRUCTION in_B={op_copy,data_stack.Last(),tos,-1};
            B->Append(in_B);
            op_stack.Append(OP_PARSE_INFO(op_nop,prec_cond,0));
            continue;}

        // Number
        if(isdigit(*str) || *str=='.'){
            char* endptr=0;
            double d=strtod(str,&endptr);
            if(!num_ok){
                LOG::cout<<"Number unexpected"<<std::endl;
                PHYSBAM_FATAL_ERROR("");}
            str=endptr;
            data_stack.Append(Add_Constant(d));
            num_ok=false;
            continue;}

        // Variable or operator
        if(isalpha(*str) || *str=='_'){
            int len=strspn(str, "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_");
            std::string var(str,len);
            str+=len;

            if(Try_Operator(block_stack,data_stack,op_stack,num_ok,var))
                continue;

            op_type type;
            if(op_tokens.func_lookup.Get(var,type)){
                if(!num_ok){
                    LOG::cout<<"Function unexpected"<<std::endl;
                    PHYSBAM_FATAL_ERROR();}
                Unwind_To_Prec(block_stack,data_stack,op_stack,prec_func);
                num_ok=false;

                OP_PARSE_INFO opi={type,prec_func,0};
                op_stack.Append(opi);
                continue;}

            int reg;
            if(!dict.Get(var,reg)){
                reg=num_tmp++;
                dict.Set(var,reg);}
            data_stack.Append(reg);
            num_ok=false;
            continue;}

        // Simple operators
        if(strchr(op_chars,*str)){
            if(str[1] && strchr(op_chars,str[1])){
                char var[3]={*str,str[1],0};
                if(Try_Operator(block_stack,data_stack,op_stack,num_ok,var)){
                    str+=2;
                    continue;}}
            char var[3]={*str,0};
            if(Try_Operator(block_stack,data_stack,op_stack,num_ok,var)){
                str++;
                continue;}}

        LOG::printf("Unrecognized character found: '%c'   rest: %s\n",*str,str);
        PHYSBAM_FATAL_ERROR();}

    PHYSBAM_ASSERT(!op_stack.m);
    PHYSBAM_ASSERT(!data_stack.m);
    PHYSBAM_ASSERT(block_stack.m==1);
    last_block=block_stack.Pop()->id;
}
const char* messages[op_last]={
    "nop\n",
    "copy  %c%d, %c%d\n",
    "add   %c%d, %c%d, %c%d\n",
    "sub   %c%d, %c%d, %c%d\n",
    "mul   %c%d, %c%d, %c%d\n",
    "div   %c%d, %c%d, %c%d\n",
    "neg   %c%d, %c%d\n",
    "inv   %c%d, %c%d\n",
    "sqrt  %c%d, %c%d\n",
    "exp   %c%d, %c%d\n",
    "ln    %c%d, %c%d\n",
    "pow   %c%d, %c%d, %c%d\n",
    "lt    %c%d, %c%d, %c%d\n",
    "le    %c%d, %c%d, %c%d\n",
    "gt    %c%d, %c%d, %c%d\n",
    "ge    %c%d, %c%d, %c%d\n",
    "eq    %c%d, %c%d, %c%d\n",
    "ne    %c%d, %c%d, %c%d\n",
    "not   %c%d, %c%d\n",
    "or    %c%d, %c%d, %c%d\n",
    "and   %c%d, %c%d, %c%d\n",
    "br_z  L%.0s%d, %c%d\n",
    "br_nz L%.0s%d, %c%d\n",
    "jmp   L%.0s%d\n",
    "label L%.0s%d\n",
    "phi   %c%d, %c%d, %c%d\n"
};
//#####################################################################
// Function Print
//#####################################################################
template<class T> void PROGRAM<T>::
Print(const INSTRUCTION& I) const
{
    printf(messages[I.type],"rcio"[I.dest>>mem_shift],I.dest&~mem_mask,
        "rcio"[I.src0>>mem_shift],I.src0&~mem_mask,
        "rcio"[I.src1>>mem_shift],I.src1&~mem_mask);
}
//#####################################################################
// Function Print
//#####################################################################
template<class T> void PROGRAM<T>::
Print(const ARRAY<INSTRUCTION>& code) const
{
    for(int i=0;i<code.m;i++){
        printf("% 3d  ",i);
        Print(code(i));}
}
//#####################################################################
// Function Print
//#####################################################################
template<class T> void PROGRAM<T>::
Print(const CODE_BLOCK* B) const
{
    for(CODE_BLOCK_NODE* N=B->head;N;N=N->next)
        Print(N->inst);
}
//#####################################################################
// Function Print
//#####################################################################
template<class T> void PROGRAM<T>::
Print() const
{
    if(!debug) return;
    for(int i=0;i<constants.m;i++)
        printf("const c%d = %g\n", i, constants(i));

    printf("block code:\n");
    for(int i=0;i<code_blocks.m;i++){
        printf("== block %i ==\n",i);
        if(!code_blocks(i)) continue;
        if(code_blocks(i)->prev[0])
            printf(messages[op_label],'L',i);
        Print(code_blocks(i));
        if(code_blocks(i)->next[0])
            printf(messages[op_jmp],'L',code_blocks(i)->next[0]->id);}

    printf("\nflat code:\n");
    Print(flat_code);
    printf("\n");
}
//#####################################################################
// Function Finalize
//#####################################################################
template<class T> void PROGRAM<T>::
Finalize()
{
    for(int j=0;j<defs(3).m;j++){
        if(!defs(3)(j)) continue;
        INSTRUCTION& o=defs(3)(j)->inst;
        if(o.type==op_copy && (o.src0&mem_mask)==mem_reg){
            INSTRUCTION& o2=Get_Def(o.src0)->inst;
            Change_Instruction(Get_Def(o.src0),o2.type,o.dest,o2.src0,o2.src1);
            Propagate_Copy(o.src0,o.dest);
            Change_Instruction(defs(3)(j),op_nop,-1,-1,-1);}}

    Leave_SSA();

    int off_in=num_tmp,off_out=off_in+var_in.m,off_const=off_out+var_out.m+extra_out;
    ARRAY<int> labels(code_blocks.m+1);
    ARRAY<CODE_BLOCK*> worklist;
    worklist.Append(code_blocks(0));
    flat_code.Remove_All();
    ARRAY<bool> done(code_blocks.m);
    done(0)=true;
    while(worklist.m){
        CODE_BLOCK* B=worklist.Pop();
        printf("%p -> %i\n",B,B->id);
        labels(B->id)=flat_code.m;
        for(CODE_BLOCK_NODE* N=B->head;N;N=N->next){
            INSTRUCTION o=N->inst;
            if(o.type==op_nop) continue;
            if(op_flags_table[o.type]&flag_reg_dest){
                if((o.dest&mem_mask)==mem_in) o.dest=(o.dest&~mem_mask)+off_in;
                else if((o.dest&mem_mask)==mem_out) o.dest=(o.dest&~mem_mask)+off_out;
                else if((o.dest&mem_mask)==mem_const) o.dest=(o.dest&~mem_mask)+off_const;}
            if(op_flags_table[o.type]&flag_reg_src0){
                if((o.src0&mem_mask)==mem_in) o.src0=(o.src0&~mem_mask)+off_in;
                else if((o.src0&mem_mask)==mem_out) o.src0=(o.src0&~mem_mask)+off_out;
                else if((o.src0&mem_mask)==mem_const) o.src0=(o.src0&~mem_mask)+off_const;}
            if(op_flags_table[o.type]&flag_reg_src1){
                if((o.src1&mem_mask)==mem_in) o.src1=(o.src1&~mem_mask)+off_in;
                else if((o.src1&mem_mask)==mem_out) o.src1=(o.src1&~mem_mask)+off_out;
                else if((o.src1&mem_mask)==mem_const) o.src1=(o.src1&~mem_mask)+off_const;}
            flat_code.Append(o);}
        if(B->next[1] && !done(B->next[1]->id)){
            done(B->next[1]->id)=true;
            worklist.Append(B->next[1]);}
        if(!B->next[0]){
            INSTRUCTION in={op_jmp,code_blocks.m,0,0};
            flat_code.Append(in);}
        else if(done(B->next[0]->id)){
            INSTRUCTION in={op_jmp,B->next[0]->id,0,0};
            flat_code.Append(in);}
        else{
            done(B->next[0]->id)=true;
            worklist.Append(B->next[0]);}}
    if(flat_code.m && flat_code.Last().type==op_jmp && flat_code.Last().dest==code_blocks.m)
        flat_code.Pop();
    labels.Last()=flat_code.m;
    for(int ip=0;ip<flat_code.m;ip++){
        INSTRUCTION& o=flat_code(ip);
        if(o.type==op_br_z || o.type==op_br_nz || o.type==op_jmp)
            o.dest=labels(o.dest);}
    finalized=true;
}
//#####################################################################
// Function Optimize
//#####################################################################
template<class T> void PROGRAM<T>::
Optimize()
{
    Sparse_Conditional_Constant_Propagation();
    Simplify_Phis();
    Copy_Propagation();
    Simplify_With_Distributive_Law();
    Local_Common_Expresssion_Elimination();
    Remove_Dead_Code();
    Compress_Registers();
}
//#####################################################################
// Function Optimize
//#####################################################################
template<class T> void PROGRAM<T>::
Make_SSA()
{
    ARRAY<ARRAY<int> > graph_in(code_blocks.m),graph_out(code_blocks.m);
    for(int i=0;i<code_blocks.m;i++){
        if(!code_blocks(i)) continue;
        for(int s=0;s<2;s++){
            if(code_blocks(i)->prev[s])
                graph_in(i).Append(code_blocks(i)->prev[s]->id);
            if(code_blocks(i)->next[s])
                graph_out(i).Append(code_blocks(i)->next[s]->id);}}

    DOMINATORS dom(graph_in,graph_out);
    dom.Compute_Frontier();

    if(debug) LOG::cout<<dom.idom<<std::endl;

    if(debug) LOG::cout<<graph_out<<std::endl;

    if(debug) LOG::cout<<dom.frontier<<std::endl;

    ARRAY<int> has_already(code_blocks.m),work(code_blocks.m),V;
    has_already.Fill(-1);
    work.Fill(-1);
    HASHTABLE<int,HASHTABLE<int> > HA;
    V=IDENTITY_ARRAY<>(num_tmp);
    HASHTABLE<int> V_used;
    V_used.Set_All(V);
    for(int bl=0;bl<code_blocks.m;bl++)
        if(CODE_BLOCK* B=code_blocks(bl))
            for(CODE_BLOCK_NODE* N=B->head;N;N=N->next){
                const INSTRUCTION& o=N->inst;
                if(op_flags_table[o.type]&flag_reg_dest)
                    HA.Get_Or_Insert(o.dest).Set(bl);
                if(op_flags_table[o.type]&flag_reg_src0)
                    HA.Get_Or_Insert(o.src0).Set(bl);
                if(op_flags_table[o.type]&flag_reg_src1)
                    HA.Get_Or_Insert(o.src1).Set(bl);}

    for(int vi=0;vi<V.m;vi++){
        int v=V(vi);
        ARRAY<int> W;
        HA.Get_Or_Insert(v).Get_Keys(W);
        work.Subset(W).Fill(vi);
        while(W.m){
            int x=W.Pop();
            for(int i=0;i<dom.frontier(x).m;i++){
                int y=dom.frontier(x)(i);
                if(has_already(y)<vi){
                    if(code_blocks(y)->prev[0] && code_blocks(y)->prev[1]){
                        INSTRUCTION in={op_phi,v,v,v};
                        code_blocks(y)->Prepend(in);}
                    has_already(y)=vi;
                    if(work(y)<vi){
                        work(y)=vi;
                        W.Append(y);}}}}}

    ARRAY<ARRAY<int> > S(V.m);
    int undef=Add_Constant(0);
    for(int i=0;i<S.m;i++)
        S(i).Append(undef);
    Make_SSA_Relabel(num_tmp,S,V_used,dom.children,0);
    for(int i=0;i<S.m;i++){
        assert(S(i).m==1 && S(i)(0)==undef);}

    if(debug) LOG::cout<<__FUNCTION__<<std::endl;
    Update_Use_Def();
    Print();
}
//#####################################################################
// Function Make_SSA_Relabel
//#####################################################################
template<class T> void PROGRAM<T>::
Make_SSA_Relabel(int& count,ARRAY<ARRAY<int> >& S,HASHTABLE<int>& V_used,const ARRAY<ARRAY<int> >& children,int bl)
{
    if(debug) LOG::cout<<count<<std::endl;
    ARRAY<int> S_pop;
    CODE_BLOCK* B=code_blocks(bl);
    for(CODE_BLOCK_NODE* N=B->head;N;N=N->next){
        INSTRUCTION& o=N->inst;
        if(op_flags_table[o.type]&flag_reg_dest){
            if(V_used.Contains(o.dest)){
                int n=count++;
                S_pop.Append(o.dest);
                S(o.dest).Append(n);
                o.dest=n;}}
        if(o.type!=op_phi){
            if(op_flags_table[o.type]&flag_reg_src0)
                if(V_used.Contains(o.src0))
                    o.src0=S(o.src0).Last();
            if(op_flags_table[o.type]&flag_reg_src1)
                if(V_used.Contains(o.src1))
                    o.src1=S(o.src1).Last();}}

    for(int s=0;s<2;s++)
        if(CODE_BLOCK* nb=code_blocks(bl)->next[s]){
            int j=nb->prev[1]==code_blocks(bl);
            for(CODE_BLOCK_NODE* N=nb->head;N && N->inst.type==op_phi;N=N->next){
                int* reg=j?&N->inst.src1:&N->inst.src0;
                *reg=S(*reg).Last();}}

    for(int i=0;i<children(bl).m;i++)
        Make_SSA_Relabel(count,S,V_used,children,children(bl)(i));

    for(int i=0;i<S_pop.m;i++)
        S(S_pop(i)).Pop();
}
//#####################################################################
// Function Leave_SSA
//#####################################################################
template<class T> void PROGRAM<T>::
Leave_SSA()
{
    ARRAY<int> var_map((IDENTITY_ARRAY<>(num_tmp)));
    // Be naive for now.
    for(int i=0;i<code_blocks.m;i++)
        if(CODE_BLOCK* B=code_blocks(i))
            for(CODE_BLOCK_NODE* N=B->head;N;N=N->next){
                INSTRUCTION& o=N->inst;
                if(o.type!=op_phi) break;
                if(o.src0==o.src1){
                    if((o.dest&mem_mask)==mem_reg) var_map(o.dest)=o.src0;
                    else{
                        o.type=op_copy;
                        continue;}}
                else{
                    INSTRUCTION in0={op_copy,o.dest,o.src0,-1};
                    INSTRUCTION in1={op_copy,o.dest,o.src1,-1};
                    code_blocks(i)->prev[0]->Append(in0);
                    code_blocks(i)->prev[1]->Append(in1);}
                o.type=op_nop;}
    Relabel_Registers(var_map);

    if(debug) LOG::cout<<__FUNCTION__<<std::endl;
    Print();
}
//#####################################################################
// Function Update_Use_Def
//#####################################################################
template<class T> void PROGRAM<T>::
Update_Use_Def()
{
    for(int i=0;i<4;i++) uses(i).Remove_All();
    uses(0).Resize(num_tmp);
    uses(1).Resize(constants.m);
    uses(2).Resize(var_in.m);
    uses(3).Resize(var_out.m+extra_out);
    defs(0)=CONSTANT_ARRAY<CODE_BLOCK_NODE*>(num_tmp,0);
    defs(3)=CONSTANT_ARRAY<CODE_BLOCK_NODE*>(var_out.m+extra_out,0);

    for(int bl=0;bl<code_blocks.m;bl++)
        if(CODE_BLOCK* B=code_blocks(bl))
            for(CODE_BLOCK_NODE* N=B->head;N;N=N->next)
                Add_Use_Def(N);
}
//#####################################################################
// Function Relabel_Registers
//#####################################################################
template<class T> void PROGRAM<T>::
Relabel_Registers(ARRAY<int>& var_map)
{
    UNION_FIND<> uf(var_map.m);
    ARRAY<int> new_map(CONSTANT_ARRAY<int>(var_map.m,-1));
    for(int i=0;i<var_map.m;i++)
        if((var_map(i)&mem_mask)==mem_reg)
            uf.Union(var_map(i),i);
    for(int i=0;i<var_map.m;i++)
        if((var_map(i)&mem_mask)!=mem_reg)
            new_map(uf.Find(i))=var_map(i);

    int k=0;
    for(int i=0;i<var_map.m;i++){
        int p=uf.Find(i);
        if(new_map(p)==-1) new_map(p)=k++;
        new_map(i)=new_map(p);}
    if(debug) LOG::cout<<new_map<<std::endl;

    for(int bl=0;bl<code_blocks.m;bl++)
        if(CODE_BLOCK* B=code_blocks(bl))
            for(CODE_BLOCK_NODE* N=B->head;N;N=N->next){
                INSTRUCTION& o=N->inst;
                Remove_Use_Def(N);
                if(op_flags_table[o.type]&flag_reg_dest && (o.dest&mem_mask)==mem_reg)
                    o.dest=new_map(o.dest);
                if(op_flags_table[o.type]&flag_reg_src0 && (o.src0&mem_mask)==mem_reg)
                    o.src0=new_map(o.src0);
                if(op_flags_table[o.type]&flag_reg_src1 && (o.src1&mem_mask)==mem_reg)
                    o.src1=new_map(o.src1);
                Add_Use_Def(N);}

    num_tmp=k;

    var_map.Exchange(new_map);
    Update_Use_Def();
}
//#####################################################################
// Function SCCP_Visit_Instruction
//#####################################################################
template<class T> void PROGRAM<T>::
SCCP_Visit_Instruction(CODE_BLOCK_NODE* N,VECTOR<ARRAY<VARIABLE_STATE>,4>& variable_state,
    ARRAY<CODE_BLOCK*>& block_worklist,ARRAY<CODE_BLOCK_NODE*>& op_worklist,ARRAY<bool>& block_exec)
{
    INSTRUCTION& o=N->inst;
    printf("try : ");Print(o);
    if(o.type==op_phi){
        VARIABLE_STATE& d=variable_state(o.dest>>mem_shift)(o.dest&~mem_mask);
        if(d.state==VARIABLE_STATE::overdetermined) return;
        VARIABLE_STATE s0=variable_state(o.src0>>mem_shift)(o.src0&~mem_mask);
        VARIABLE_STATE s1=variable_state(o.src1>>mem_shift)(o.src1&~mem_mask);

        printf("sccp %i %g : %i %g %i %g : ",d.state,d.value,s0.state,s0.value,s1.state,s1.value);Print(o);

        if(s0.state==VARIABLE_STATE::overdetermined || s1.state==VARIABLE_STATE::overdetermined){
            d.state=VARIABLE_STATE::overdetermined;
            op_worklist.Append(N);
            return;}

        if(s0.state==VARIABLE_STATE::constant && s1.state==VARIABLE_STATE::constant){
            if(s0.value!=s1.value) d.state=VARIABLE_STATE::overdetermined;
            else if(d.state==VARIABLE_STATE::constant) return;
            else{d.state=VARIABLE_STATE::constant;d.value=s0.value;}
            op_worklist.Append(N);
            return;}

        if(s0.state==VARIABLE_STATE::constant || s1.state==VARIABLE_STATE::constant){
            if(d.state==VARIABLE_STATE::constant) return;
            d.state=VARIABLE_STATE::constant;
            d.value=s0.state==VARIABLE_STATE::constant?s0.value:s1.value;
            op_worklist.Append(N);
            return;}
        return;}

    if(op_flags_table[o.type]&flag_reg_dest){
        VARIABLE_STATE& d=variable_state(o.dest>>mem_shift)(o.dest&~mem_mask);
        if(d.state==VARIABLE_STATE::overdetermined) return;
        VARIABLE_STATE* s0=op_flags_table[o.type]&flag_reg_src0?&variable_state(o.src0>>mem_shift)(o.src0&~mem_mask):0;
        VARIABLE_STATE* s1=op_flags_table[o.type]&flag_reg_src1?&variable_state(o.src1>>mem_shift)(o.src1&~mem_mask):0;
        printf("sccp %i %g : %p %p : ",d.state,d.value,s0,s1);Print(o);

        // rhs is known; evaluate instruction
        if((!s0 || s0->state==VARIABLE_STATE::constant) && (!s1 || s1->state==VARIABLE_STATE::constant)){
            if(d.state==VARIABLE_STATE::constant) return;
            printf("constant ");Print(o);
            d.state=VARIABLE_STATE::constant;
            d.value=Evaluate_Op(o.type,s0?s0->value:0,s1?s1->value:0);
            op_worklist.Append(N);
            return;}

        // Might be able to deduce result with partial information
        if(op_flags_table[o.type]&flag_can_deduce){
            // rhs is partially known; see if result can be deduced anyway
            if((s0 && s0->state==VARIABLE_STATE::constant) || (s1 && s1->state==VARIABLE_STATE::constant)){
                T* a=(s0 && s0->state==VARIABLE_STATE::constant)?&s0->value:0;
                T* b=(s1 && s1->state==VARIABLE_STATE::constant)?&s1->value:0;
                if(Deduce_Op(o.type,d.value,a,b)){
                    if(d.state==VARIABLE_STATE::constant) return;
                    printf("constant ");Print(o);
                    d.state=VARIABLE_STATE::constant;
                    op_worklist.Append(N);
                    return;}
                else if((s0 && s0->state==VARIABLE_STATE::overdetermined) || (s1 && s1->state==VARIABLE_STATE::overdetermined)){
                    printf("overdetermined ");Print(o);
                    d.state=VARIABLE_STATE::overdetermined;
                    op_worklist.Append(N);
                    return;}}
            if((!s0 || s0->state==VARIABLE_STATE::overdetermined) && (!s1 || s1->state==VARIABLE_STATE::overdetermined)){
                printf("overdetermined ");Print(o);
                d.state=VARIABLE_STATE::overdetermined;
                op_worklist.Append(N);
                return;}}
        else if((s0 && s0->state==VARIABLE_STATE::overdetermined) || (s1 && s1->state==VARIABLE_STATE::overdetermined)){
            printf("overdetermined ");Print(o);
            d.state=VARIABLE_STATE::overdetermined;
            op_worklist.Append(N);
            return;}
        return;}

    if(o.type==op_jmp){
        if(!block_exec(o.dest)){
            block_exec(o.dest)=true;
            block_worklist.Append(code_blocks(o.dest));}
        return;}

    if(o.type==op_nop || o.type==op_label) return;
    // This leaves only op_br_z and op_br_nz

    VARIABLE_STATE s0=variable_state(o.src0>>mem_shift)(o.src0&~mem_mask);
    if(s0.state==VARIABLE_STATE::overdetermined){
        if(!block_exec(o.dest)){
            block_exec(o.dest)=true;
            block_worklist.Append(code_blocks(o.dest));}
        if(!block_exec(o.src1)){
            block_exec(o.src1)=true;
            block_worklist.Append(code_blocks(o.src1));}
        return;}

    if(s0.state!=VARIABLE_STATE::constant) return;

    int dst=o.dest;
    if((o.type==op_br_nz)!=(s0.value!=0)) dst=o.src1;
    if(!block_exec(dst)){
        block_exec(dst)=true;
        block_worklist.Append(code_blocks(dst));}
}
//#####################################################################
// Function Sparse_Conditional_Constant_Propagation
//#####################################################################
template<class T> void PROGRAM<T>::
Sparse_Conditional_Constant_Propagation()
{
    VECTOR<ARRAY<VARIABLE_STATE>,4> variable_state;

    VARIABLE_STATE unknown={0,VARIABLE_STATE::unknown};
    VARIABLE_STATE overdetermined={0,VARIABLE_STATE::overdetermined};
    VARIABLE_STATE constant={0,VARIABLE_STATE::constant};
    variable_state(mem_reg>>mem_shift)=CONSTANT_ARRAY<VARIABLE_STATE>(num_tmp,unknown);
    variable_state(mem_out>>mem_shift)=CONSTANT_ARRAY<VARIABLE_STATE>(var_out.m+extra_out,unknown);
    variable_state(mem_in>>mem_shift)=CONSTANT_ARRAY<VARIABLE_STATE>(var_in.m,overdetermined);
    variable_state(mem_const>>mem_shift)=CONSTANT_ARRAY<VARIABLE_STATE>(constants.m,constant);

    for(int i=0;i<constants.m;i++)
        variable_state(mem_const>>mem_shift)(i).value=constants(i);

    ARRAY<CODE_BLOCK*> block_worklist;
    ARRAY<CODE_BLOCK_NODE*> op_worklist;
    ARRAY<bool> block_exec(code_blocks.m);
    block_exec(0)=true;
    block_worklist.Append(code_blocks(0));

    while(block_worklist.m || op_worklist.m){
        while(op_worklist.m){
            CODE_BLOCK_NODE* N=op_worklist.Pop();
            for(typename HASHTABLE<CODE_BLOCK_NODE*>::ITERATOR it(Get_Uses(N));it.Valid();it.Next()){
                if(block_exec(it.Key()->block->id))
                    SCCP_Visit_Instruction(it.Key(),variable_state,block_worklist,op_worklist,block_exec);}}

        while(block_worklist.m){
            CODE_BLOCK* block=block_worklist.Pop();
            for(CODE_BLOCK_NODE* N=block->head;N;N=N->next)
                SCCP_Visit_Instruction(N,variable_state,block_worklist,op_worklist,block_exec);
            if(block->next[0] && (!block->tail || (block->tail->inst.type!=op_jmp && block->tail->inst.type!=op_br_nz && block->tail->inst.type!=op_br_z)))
                if(!block_exec(block->next[0]->id)){
                    block_exec(block->next[0]->id)=true;
                    block_worklist.Append(block->next[0]);}}}

    for(int i=0;i<4;i++)
        for(int j=0;j<defs(i).m;j++)
            if(variable_state(i)(j).state==VARIABLE_STATE::constant)
                if(CODE_BLOCK_NODE* N=defs(i)(j)){
                    INSTRUCTION& o=N->inst;
                    if(debug) LOG::cout<<"solved:   r"<<j<<" = "<<variable_state(i)(j).value<<std::endl;
                    Change_Instruction(N,op_copy,o.dest,Add_Constant(variable_state(i)(j).value),-1);
                    Propagate_Copy_And_Remove_Instruction(N);}

    for(int i=0;i<code_blocks.m;i++)
        if(CODE_BLOCK* B=code_blocks(i))
            if(B->tail){
                INSTRUCTION& o=B->tail->inst;
                if((o.type==op_br_z || o.type==op_br_z) && (o.src0&mem_mask)==mem_const){
                    printf("eliminate %i %g ",o.src0&~mem_mask,constants(o.src0&~mem_mask));Print(o);
                    if((o.type==op_br_nz)==(constants(o.src0&~mem_mask)!=0)){
                        o.dest=o.src1;
                        exchange(B->next[0],B->next[1]);}
                    Change_Instruction(B->tail,op_nop,-1,-1,-1);
                    Detach_Next(B,1);}}

    if(debug) LOG::cout<<__FUNCTION__<<std::endl;
    Print();
}
//#####################################################################
// Function Copy_Propagation
//#####################################################################
template<class T> void PROGRAM<T>::
Propagate_Copy(int old_var,int new_var)
{
    assert(old_var!=new_var);
    HASHTABLE<CODE_BLOCK_NODE*>& ls=Get_Uses(old_var);
    while(ls.Size()){
        typename HASHTABLE<CODE_BLOCK_NODE*>::ITERATOR it(ls);
        CODE_BLOCK_NODE* N=it.Key();
        int flags=op_flags_table[N->inst.type];
        bool d0=flags&flag_reg_src0 && N->inst.src0==old_var;
        bool d1=flags&flag_reg_src1 && N->inst.src1==old_var;
        assert(d0 || d1);
        Remove_Use_Def(N);
        if(d0) N->inst.src0=new_var;
        if(d1) N->inst.src1=new_var;
        Add_Use_Def(N);}
}
//#####################################################################
// Function Copy_Propagation
//#####################################################################
template<class T> void PROGRAM<T>::
Copy_Propagation()
{
    // Copy propagation
    for(int bl=0;bl<code_blocks.m;bl++)
        if(CODE_BLOCK* B=code_blocks(bl))
            for(CODE_BLOCK_NODE* N=B->head;N;){
                CODE_BLOCK_NODE* N2=N;
                N=N->next;
                if(N2->inst.type==op_copy || (N2->inst.type==op_phi && N2->inst.src0==N2->inst.src1))
                    Propagate_Copy_And_Remove_Instruction(N2);}

    if(debug) LOG::cout<<__FUNCTION__<<std::endl;
    Print();
}
//#####################################################################
// Function Copy_Propagation
//#####################################################################
// Don't assume o.type==op_copy, but do assume it behaves as dest=src0
template<class T> void PROGRAM<T>::
Propagate_Copy_And_Remove_Instruction(CODE_BLOCK_NODE* N)
{
    INSTRUCTION& o=N->inst;
    Propagate_Copy(o.dest,o.src0);
    if((N->inst.dest&mem_mask)==mem_reg) Delete_Instruction(N);
    else Change_Instruction(N,op_copy,o.dest,o.src0,-1);
}
//#####################################################################
// Function Simplify_Phis
//#####################################################################
template<class T> void PROGRAM<T>::
Simplify_Phis()
{
    if(debug) LOG::cout<<__FUNCTION__<<std::endl;
    Print();
    // Replace phi's with unreachble input by a copy.
    for(int bl=0;bl<code_blocks.m;bl++)
        if(CODE_BLOCK* B=code_blocks(bl)){
            if(!B->prev[0] || !B->prev[1]){
                int j=!B->prev[0];
                for(CODE_BLOCK_NODE* N=B->head;N;){
                    INSTRUCTION& o=N->inst;
                    if(o.type!=op_phi) break;
                    CODE_BLOCK_NODE* N2=N;
                    N=N->next;
                    Change_Instruction(N2,op_copy,o.dest,j?o.src1:o.src0,-1);
                    Propagate_Copy_And_Remove_Instruction(N2);}
                B->prev[0]=B->prev[j];}
            else{
                for(CODE_BLOCK_NODE* N=B->head;N;){
                    INSTRUCTION& o=N->inst;
                    if(o.type!=op_phi) break;
                    CODE_BLOCK_NODE* N2=N;
                    N=N->next;
                    if(o.src0!=o.src1) continue;
                    Propagate_Copy_And_Remove_Instruction(N2);}}}
    
    if(debug) LOG::cout<<__FUNCTION__<<std::endl;
    Print();
    if(debug) LOG::cout<<__FUNCTION__<<std::endl;
}
//#####################################################################
// Function Remove_Dead_Code
//#####################################################################
template<class T> void PROGRAM<T>::
Remove_Dead_Code()
{
    // Prune unreachable basic blocks (bl=0 has no prev)
    ARRAY<CODE_BLOCK*> worklist;
    for(int bl=1;bl<code_blocks.m;bl++)
        if(CODE_BLOCK* B=code_blocks(bl))
            if(!B->prev[0] && !B->prev[1])
                worklist.Append(B);
    for(int i=0;i<worklist.m;i++){
        CODE_BLOCK* B=worklist(i);
        for(int s=0;s<2;s++)
            if(CODE_BLOCK* C=B->next[s]){
                Detach_Next(B,s);
                if(!C->prev[0] && !C->prev[1])
                    worklist.Append(C);}}
    Simplify_Phis();
    for(int i=0;i<worklist.m;i++){
        CODE_BLOCK* B=worklist(i);
        while(B->head) Delete_Instruction(B->head);
        code_blocks(B->id)=0;
        delete B;}

    // Merge basic blocks
    for(int bl=0;bl<code_blocks.m;bl++)
        if(CODE_BLOCK* B=code_blocks(bl))
            while(B->next[0] && !B->next[1] && !B->next[0]->prev[1]){
                CODE_BLOCK* nb=B->next[0];
                for(int s=0;s<2;s++){
                    if(nb->next[s]){
                        int j=nb->next[s]->prev[1]==nb;
                        nb->next[s]->prev[j]=B;}
                    B->next[s]=nb->next[s];}
                B->Merge_Code(nb);
                delete code_blocks(nb->id);
                code_blocks(nb->id)=0;}

    Update_Use_Def();
    ARRAY<CODE_BLOCK_NODE*> worklist2;
    for(int j=0;j<defs(0).m;j++)
        if(CODE_BLOCK_NODE* N=defs(0)(j))
            if(!uses(0)(j).Size())
                worklist2.Append(N);

    while(worklist2.m){
        CODE_BLOCK_NODE* N=worklist2.Pop();
        int flags=op_flags_table[N->inst.type],r0=N->inst.src0,r1=N->inst.src1;
        if(flags&flag_reg_src0 && (r0&mem_mask)==mem_reg && Get_Uses(r0).Size()==1) worklist2.Append(Get_Def(r0));
        if(flags&flag_reg_src1 && (r1&mem_mask)==mem_reg && r1!=r0 && Get_Uses(r1).Size()==1) worklist2.Append(Get_Def(r1));
        Delete_Instruction(N);}

    for(int bl=0;bl<code_blocks.m;bl++)
        if(CODE_BLOCK* B=code_blocks(bl))
            for(CODE_BLOCK_NODE* N=B->head;N;){
                CODE_BLOCK_NODE* N2=N;
                N=N->next;
                if(N2->inst.type==op_label || N2->inst.type==op_nop)
                    Delete_Instruction(N2);}

    if(debug) LOG::cout<<__FUNCTION__<<std::endl;
    Print();
}
//#####################################################################
// Function Detach_Next
//#####################################################################
template<class T> void PROGRAM<T>::
Detach_Next(CODE_BLOCK* B,int j)
{
    int phi_j=B->next[j]->prev[1]==B;
    B->next[j]->prev[phi_j]=0;
    B->next[j]=0;
}
//#####################################################################
// Function Add_Constant
//#####################################################################
template<class T> int PROGRAM<T>::
Add_Constant(T x)
{
    typename std::map<T,int>::iterator it=constant_lookup.find(x);
    if(it!=constant_lookup.end()) return it->second|mem_const;
    int i=constants.Append(x);
    uses(1).Append(HASHTABLE<CODE_BLOCK_NODE*>());
    constant_lookup[x]=i;
    return i|mem_const;
}
//#####################################################################
// Function Local_Common_Expresssion_Elimination
//#####################################################################
template<class T> void PROGRAM<T>::
Local_Common_Expresssion_Elimination(CODE_BLOCK* B)
{
    HASHTABLE<VECTOR<int,3>,int> expr;
    for(CODE_BLOCK_NODE* N=B->head;N;N=N->next){
        INSTRUCTION& o=N->inst;
        if(!(op_flags_table[o.type]&flag_reg_dest)) continue;
        Reduce_In_Place(N);
        if(o.type==op_copy){
            Propagate_Copy(o.dest,o.src0);
            continue;}
        VECTOR<int,3> in(o.type,o.src0,o.src1);
        if(!(op_flags_table[o.type]&flag_reg_src1)) in.z=-1;
        if(op_flags_table[o.type]&flag_commute) exchange_sort(in.y,in.z);
        int e=-1;
        if(expr.Get(in,e))
            Propagate_Copy(o.dest,e);
        else expr.Set(in,o.dest);}
}
//#####################################################################
// Function Local_Common_Expresssion_Elimination
//#####################################################################
template<class T> void PROGRAM<T>::
Local_Common_Expresssion_Elimination()
{
    if(debug) LOG::cout<<__FUNCTION__<<std::endl;
    for(int bl=0;bl<code_blocks.m;bl++){
        if(!code_blocks(bl)) continue;
        Local_Common_Expresssion_Elimination(code_blocks(bl));}
    Print();
    if(debug) LOG::cout<<__FUNCTION__<<std::endl;
}
//#####################################################################
// Function Reduce_In_Place
//#####################################################################
template<class T> void PROGRAM<T>::
Reduce_In_Place(CODE_BLOCK_NODE* N)
{
    int const_0=Add_Constant(0),const_1=Add_Constant(1);
    INSTRUCTION& o=N->inst;
    switch(o.type){
        case op_add:
            if(o.src1==const_0){Change_Instruction(N,op_copy,o.dest,o.src0,-1);break;}
            if(o.src0==const_0){Change_Instruction(N,op_copy,o.dest,o.src1,-1);break;}
            break;
        case op_sub:
            if(o.src0==o.src1){Change_Instruction(N,op_copy,o.dest,const_0,-1);break;}
            if(o.src1==const_0){Change_Instruction(N,op_copy,o.dest,o.src0,-1);break;}
            if(o.src0==const_0){Change_Instruction(N,op_neg,o.dest,o.src1,-1);break;}
            break;
        case op_mul:
            if(o.src1==const_0 || o.src0==const_0){Change_Instruction(N,op_copy,o.dest,const_0,-1);break;}
            if(o.src1==const_1){Change_Instruction(N,op_copy,o.dest,o.src0,-1);break;}
            if(o.src0==const_1){Change_Instruction(N,op_copy,o.dest,o.src1,-1);break;}
            if(o.src1==Add_Constant(2)){Change_Instruction(N,op_add,o.dest,o.src0,o.src0);break;}
            if(o.src0==Add_Constant(2)){Change_Instruction(N,op_add,o.dest,o.src1,o.src1);break;}
            break;
        case op_div:
            if(o.src0==o.src1){Change_Instruction(N,op_copy,o.dest,const_1,-1);break;}
            if(o.src0==const_0 || o.src1==const_1){Change_Instruction(N,op_copy,o.dest,o.src0,-1);break;}
            if(o.src0==const_1){Change_Instruction(N,op_inv,o.dest,o.src1,-1);break;}
            break;
        case op_neg: break;
        case op_inv: break;
        case op_sqrt: break;
        case op_exp: break;
        case op_ln: break;
        case op_pow:
            if(o.src1==const_0){Change_Instruction(N,op_copy,o.dest,const_1,-1);break;}
            if(o.src1==const_1){Change_Instruction(N,op_copy,o.dest,o.src0,-1);break;}
            if(o.src1==Add_Constant(2)){Change_Instruction(N,op_mul,o.dest,o.src0,o.src0);break;}
            if(o.src1==Add_Constant(-1)){Change_Instruction(N,op_inv,o.dest,o.src0,-1);break;}
            break;
        case op_lt:
        case op_gt:
        case op_ne:
            if(o.src0==o.src1){Change_Instruction(N,op_copy,o.dest,const_0,-1);break;}
            break;
        case op_le:
        case op_ge:
        case op_eq:
            if(o.src0==o.src1){Change_Instruction(N,op_copy,o.dest,const_1,-1);break;}
            break;
        case op_or:
        case op_and:
            if(o.src0==o.src1){Change_Instruction(N,op_copy,o.dest,o.src0,-1);break;}
            break;
        default: break;}
}
//#####################################################################
// Function Simplify_With_Distributive_Law_Helper_Prod
//#####################################################################
template<class T> void PROGRAM<T>::
Simplify_With_Distributive_Law_Helper_Prod(int var,ARRAY<CODE_BLOCK_NODE*>& extra_nodes,DISTRIBUTE_HELPER<T>& prod)
{
    if((var&mem_mask)==mem_const){
        prod.coeff*=constants(var&~mem_mask);
        return;}

    if((var&mem_mask)==mem_in || (var&mem_mask)==mem_out || Get_Uses(var).Size()>1){
        prod.product.Append(var);
        return;}

    CODE_BLOCK_NODE* N=Get_Def(var);
    assert(N);

    if(N->inst.type==op_neg){
        extra_nodes.Append(N);
        prod.coeff=-prod.coeff;
        return Simplify_With_Distributive_Law_Helper_Prod(N->inst.src0,extra_nodes,prod);}

    if(N->inst.type!=op_mul){
        prod.product.Append(var);
        return;}

    extra_nodes.Append(N);
    Simplify_With_Distributive_Law_Helper_Prod(N->inst.src0,extra_nodes,prod);
    Simplify_With_Distributive_Law_Helper_Prod(N->inst.src1,extra_nodes,prod);
}
//#####################################################################
// Function Simplify_With_Distributive_Law_Helper
//#####################################################################
template<class T> void PROGRAM<T>::
Simplify_With_Distributive_Law_Helper(int var,T coeff,ARRAY<CODE_BLOCK_NODE*>& extra_nodes,ARRAY<DISTRIBUTE_HELPER<T> >& summands)
{
    if((var&mem_mask)==mem_const){
        DISTRIBUTE_HELPER<T> prod={coeff*constants(var&~mem_mask)};
        summands.Append(prod);
        return;}

    if(extra_nodes.m && ((var&mem_mask)==mem_in || (var&mem_mask)==mem_out || Get_Uses(var).Size()>1)){
        DISTRIBUTE_HELPER<T> prod={coeff};
        summands(summands.Append(prod)).product.Append(var);
        return;}

    CODE_BLOCK_NODE* N=Get_Def(var);
    assert(N);

    if(N->inst.type==op_mul){
        if((N->inst.src0&mem_mask)==mem_const){
            extra_nodes.Append(N);
            return Simplify_With_Distributive_Law_Helper(N->inst.src1,coeff*constants(N->inst.src0&~mem_mask),extra_nodes,summands);}

        if((N->inst.src1&mem_mask)==mem_const){
            extra_nodes.Append(N);
            return Simplify_With_Distributive_Law_Helper(N->inst.src0,coeff*constants(N->inst.src1&~mem_mask),extra_nodes,summands);}

        DISTRIBUTE_HELPER<T> prod={coeff};
        return Simplify_With_Distributive_Law_Helper_Prod(var,extra_nodes,summands(summands.Append(prod)));}

    if(N->inst.type==op_neg){
        extra_nodes.Append(N);
        return Simplify_With_Distributive_Law_Helper(N->inst.src0,-coeff,extra_nodes,summands);}

    if(N->inst.type==op_copy){
        extra_nodes.Append(N);
        return Simplify_With_Distributive_Law_Helper(N->inst.src0,coeff,extra_nodes,summands);}

    if(N->inst.type!=op_add && N->inst.type!=op_sub){
        DISTRIBUTE_HELPER<T> prod={coeff};
        summands(summands.Append(prod)).product.Append(var);
        return;}

    extra_nodes.Append(N);
    Simplify_With_Distributive_Law_Helper(N->inst.src0,coeff,extra_nodes,summands);
    if(N->inst.type==op_sub) coeff=-coeff;
    Simplify_With_Distributive_Law_Helper(N->inst.src1,coeff,extra_nodes,summands);
}
//#####################################################################
// Function Simplify_With_Distributive_Law_Emit_Prod
//#####################################################################
template<class T> int PROGRAM<T>::
Simplify_With_Distributive_Law_Emit_Prod(DISTRIBUTE_HELPER<T>& prod,CODE_BLOCK_NODE* before,bool& need_neg)
{
    need_neg=false;
    if(prod.coeff==0) return Add_Constant(0);
    if(prod.product.m==0) return Add_Constant(prod.coeff);

    while(prod.product.m>1){
        int new_reg=num_tmp++,src1=prod.product.Pop();
        Add_Raw_Instruction_To_Block_Before(before->block,before,op_mul,new_reg,prod.product.Last(),src1);
        prod.product.Last()=new_reg;}

    if(prod.coeff==1) return prod.product(0);
    if(prod.coeff==-1){
        need_neg=true;
        return prod.product(0);}

    int new_reg=num_tmp++;
    Add_Raw_Instruction_To_Block_Before(before->block,before,op_mul,new_reg,Add_Constant(prod.coeff),prod.product.Last());
    return new_reg;
}
//#####################################################################
// Function Simplify_With_Distributive_Law_Emit_Simple
//#####################################################################
template<class T> int PROGRAM<T>::
Simplify_With_Distributive_Law_Emit_Simple(ARRAY<DISTRIBUTE_HELPER<T> >& summands,CODE_BLOCK_NODE* before,bool& need_neg)
{
    int accum_reg=-1;
    need_neg=false;
    while(summands.m){
        bool neg_a=false;
        int out_a=Simplify_With_Distributive_Law_Emit_Prod(summands.Last(),before,neg_a);
        summands.Pop();
        if(accum_reg==-1){
            accum_reg=out_a;
            need_neg=neg_a;}
        else{
            int new_reg=num_tmp++;
            if(need_neg==neg_a) Add_Raw_Instruction_To_Block_Before(before->block,before,op_add,new_reg,accum_reg,out_a);
            else{
                if(need_neg) exchange(accum_reg,out_a);
                Add_Raw_Instruction_To_Block_Before(before->block,before,op_sub,new_reg,accum_reg,out_a);
                need_neg=false;}
            accum_reg=new_reg;}}
    return accum_reg;
}
//#####################################################################
// Function Simplify_With_Distributive_Law_Emit_Consts
//#####################################################################
template<class T> int PROGRAM<T>::
Simplify_With_Distributive_Law_Emit_Consts(ARRAY<DISTRIBUTE_HELPER<T> >& summands,CODE_BLOCK_NODE* before,bool& need_neg)
{
    need_neg=false;
    T const_term=0;

    std::map<T,ARRAY<DISTRIBUTE_HELPER<T> > > consts;
    for(int i=0;i<summands.m;i++){
        if(summands(i).product.m==0) const_term+=summands(i).coeff;
        else consts[abs(summands(i).coeff)].Append(summands(i));}

    int accum_reg=-1;

    ARRAY<DISTRIBUTE_HELPER<T> > no_coeff;
    typename std::map<T,ARRAY<DISTRIBUTE_HELPER<T> > >::iterator it;
    if((it=consts.find(0))!=consts.end()) consts.erase(it);
    if((it=consts.find(1))!=consts.end()){
        no_coeff.Exchange(it->second);
        consts.erase(it);}

    for(auto& it:consts){
        bool neg_a=false;
        for(int i=0;i<it.second.m;i++)
            it.second(i).coeff=sign(it.second(i).coeff);
        int out_a=Simplify_With_Distributive_Law_Emit_Simple(it.second,before,neg_a);
        T c=it.first;
        if(neg_a) c=-c;
        int new_reg=num_tmp++;
        Add_Raw_Instruction_To_Block_Before(before->block,before,op_mul,new_reg,Add_Constant(c),out_a);
        if(accum_reg!=-1){
            int new_reg2=num_tmp++;
            Add_Raw_Instruction_To_Block_Before(before->block,before,op_add,new_reg2,accum_reg,new_reg);
            accum_reg=new_reg2;}
        else accum_reg=new_reg;}

    if(no_coeff.m){
        int out_a=Simplify_With_Distributive_Law_Emit_Simple(no_coeff,before,need_neg);
        if(accum_reg!=-1){
            int new_reg=num_tmp++;
            Add_Raw_Instruction_To_Block_Before(before->block,before,need_neg?op_sub:op_add,new_reg,accum_reg,out_a);
            need_neg=false;
            accum_reg=new_reg;}
        else accum_reg=out_a;}

    if(accum_reg==-1) return Add_Constant(const_term);

    if(const_term){
        int new_reg=num_tmp++;
        Add_Raw_Instruction_To_Block_Before(before->block,before,need_neg?op_sub:op_add,new_reg,Add_Constant(const_term),accum_reg);
        need_neg=false;
        accum_reg=new_reg;}
    return accum_reg;
}
//#####################################################################
// Function Simplify_With_Distributive_Law_Emit
//#####################################################################
template<class T> int PROGRAM<T>::
Simplify_With_Distributive_Law_Emit(ARRAY<DISTRIBUTE_HELPER<T> >& summands,CODE_BLOCK_NODE* before,bool& need_neg)
{
    HASHTABLE<int,HASHTABLE<int> > has_var;
    for(int i=0;i<summands.m;i++)
        for(int j=0;j<summands(i).product.m;j++)
            has_var.Get_Or_Insert(summands(i).product(j)).Set(i);

    int max_count=-1,best_var=-1;
    for(HASHTABLE<int,HASHTABLE<int> >::ITERATOR it(has_var);it.Valid();it.Next())
        if(it.Data().Size()>max_count){
            max_count=it.Data().Size();
            best_var=it.Key();}

    if(max_count<=1){
        has_var.Remove_All();
        return Simplify_With_Distributive_Law_Emit_Consts(summands,before,need_neg);}

    int k=0;
    ARRAY<DISTRIBUTE_HELPER<T> > new_summands(max_count);
    HASHTABLE<int>& h=has_var.Get(best_var);
    for(HASHTABLE<int>::ITERATOR it(h);it.Valid();it.Next()){
        int n=it.Key();
        new_summands(k).coeff=summands(n).coeff;
        new_summands(k++).product.Exchange(summands(n).product);
        summands(n).coeff=0;}
    for(int i=summands.m-1;i>=0;i--)
        if(summands(i).coeff==0){
            summands(i).coeff=summands.Last().coeff;
            summands(i).product.Exchange(summands.Last().product);
            summands.Pop();}
    for(int i=0;i<new_summands.m;i++)
        new_summands(i).product.Remove_Index_Lazy(new_summands(i).product.Find(best_var));
    has_var.Remove_All();
    bool neg_a=false,neg_b=false;
    int out_a=Simplify_With_Distributive_Law_Emit(new_summands,before,neg_a),out_reg=num_tmp++;
    Add_Raw_Instruction_To_Block_Before(before->block,before,op_mul,out_reg,best_var,out_a);
    if(!summands.m){
        need_neg=neg_a;
        return out_reg;}
    int out_b=Simplify_With_Distributive_Law_Emit(summands,before,neg_b);
    int ret_reg=num_tmp++;
    if(neg_a==neg_b){
        Add_Raw_Instruction_To_Block_Before(before->block,before,op_add,ret_reg,out_reg,out_b);
        need_neg=neg_a;
        return ret_reg;}
    need_neg=false;
    if(neg_b) Add_Raw_Instruction_To_Block_Before(before->block,before,op_sub,ret_reg,out_reg,out_b);
    else Add_Raw_Instruction_To_Block_Before(before->block,before,op_sub,ret_reg,out_b,out_reg);
    return ret_reg;
}
//#####################################################################
// Function Simplify_With_Distributive_Law
//#####################################################################
template<class T> CODE_BLOCK_NODE* PROGRAM<T>::
Simplify_With_Distributive_Law(CODE_BLOCK_NODE* N)
{
    if(debug) LOG::cout<<N<<std::endl;

    ARRAY<CODE_BLOCK_NODE*> extra_nodes;
    ARRAY<DISTRIBUTE_HELPER<T> > summands;
    Simplify_With_Distributive_Law_Helper(N->inst.dest,1,extra_nodes,summands);
    if(!extra_nodes.m) return N->prev;
    extra_nodes.Prune_Duplicates();
    for(int i=0;i<extra_nodes.m;i++)
        if(extra_nodes(i)!=N)
            Delete_Instruction(extra_nodes(i));
    bool need_neg=false;
    CODE_BLOCK_NODE* P=N->prev;
    int out_reg=Simplify_With_Distributive_Law_Emit(summands,N,need_neg);
    if(need_neg) Change_Instruction(N,op_neg,N->inst.dest,out_reg,-1);
    else{
        Change_Instruction(N,op_copy,N->inst.dest,out_reg,-1);
        Propagate_Copy_And_Remove_Instruction(N);}
    return P;
}
//#####################################################################
// Function Simplify_With_Distributive_Law
//#####################################################################
template<class T> void PROGRAM<T>::
Simplify_With_Distributive_Law()
{
    for(int bl=0;bl<code_blocks.m;bl++)
        if(CODE_BLOCK* B=code_blocks(bl))
            for(CODE_BLOCK_NODE* N=B->tail;N;){
                CODE_BLOCK_NODE* N2=N;
                N=N->prev;
                if(N2->inst.type!=op_neg && N2->inst.type!=op_add && N2->inst.type!=op_sub && N2->inst.type!=op_mul && N2->inst.type!=op_copy) continue;
                if((N2->inst.dest&mem_mask)!=mem_in && (N2->inst.dest&mem_mask)!=mem_out && Get_Uses(N2->inst.dest).Size()<=1) continue;
                N=Simplify_With_Distributive_Law(N2);}
    if(debug) LOG::cout<<__FUNCTION__<<std::endl;
    Print();
    if(debug) LOG::cout<<__FUNCTION__<<std::endl;
}
//#####################################################################
// Function Eliminate_Unused_Register
//#####################################################################
template<class T> void PROGRAM<T>::
Eliminate_Unused_Register(int k)
{
    int last=uses(0).m-1;
    if(k<last){
        if(defs(0)(last)) defs(0)(last)->inst.dest=k;
        defs(0)(k)=defs(0)(last);
        Propagate_Copy(last,k);}
    defs(0).Pop();
    uses(0).Pop();
}
//#####################################################################
// Function Eliminate_Unused_Constant
//#####################################################################
template<class T> void PROGRAM<T>::
Eliminate_Unused_Constant(int k)
{
    int last=uses(1).m-1;
    constant_lookup.erase(constant_lookup.find(constants(k)));
    if(k<last){
        constant_lookup[constants(last)]=k;
        constants(k)=constants(last);
        Propagate_Copy(last|mem_const,k|mem_const);}
    constants.Pop();
    uses(1).Pop();
}
//#####################################################################
// Function Compress_Registers
//#####################################################################
template<class T> void PROGRAM<T>::
Compress_Registers()
{
    


    for(int j=0;j<defs(0).m;j++)
        if(!defs(0)(j))
            Eliminate_Unused_Register(j--);

    for(int j=0;j<constants.m;j++)
        if(!uses(1)(j).Size())
            Eliminate_Unused_Constant(j--);

    num_tmp=uses(0).m;

    if(debug) LOG::cout<<__FUNCTION__<<std::endl;
    Print();
}
//#####################################################################
// Function Remove_Use_Def
//#####################################################################
template<class T> void PROGRAM<T>::
Remove_Use_Def(CODE_BLOCK_NODE* N)
{
    int flags=op_flags_table[N->inst.type];
    if(flags&flag_reg_dest) Set_Def(N->inst.dest,0);
    if(flags&flag_reg_src0) Get_Uses(N->inst.src0).Delete_If_Present(N);
    if(flags&flag_reg_src1) Get_Uses(N->inst.src1).Delete_If_Present(N);
}
//#####################################################################
// Function Add_Use_Def
//#####################################################################
template<class T> void PROGRAM<T>::
Add_Use_Def(CODE_BLOCK_NODE* N)
{
    int flags=op_flags_table[N->inst.type];
    if(flags&flag_reg_dest){
        ARRAY<CODE_BLOCK_NODE*>& a=defs(N->inst.dest>>mem_shift);
        int i=N->inst.dest&~mem_mask;
        if(i>=a.m) a.Resize(i+1);
        a(i)=N;}
    if(op_flags_table[N->inst.type]&flag_reg_src0){
        ARRAY<HASHTABLE<CODE_BLOCK_NODE*> >& a=uses(N->inst.src0>>mem_shift);
        int i=N->inst.src0&~mem_mask;
        if(i>=a.m) a.Resize(i+1);
        HASHTABLE<CODE_BLOCK_NODE*>& l=a(i);
        l.Set(N);}
    if(op_flags_table[N->inst.type]&flag_reg_src1){
        ARRAY<HASHTABLE<CODE_BLOCK_NODE*> >& a=uses(N->inst.src1>>mem_shift);
        int i=N->inst.src1&~mem_mask;
        if(i>=a.m) a.Resize(i+1);
        a(i).Set(N);}
}
//#####################################################################
// Function Delete_Instruction
//#####################################################################
template<class T> void PROGRAM<T>::
Delete_Instruction(CODE_BLOCK_NODE* N)
{
    Remove_Use_Def(N);
    N->block->Delete_Node(N); // Remove from list but do not delete
}
//#####################################################################
// Function Change_Instruction
//#####################################################################
template<class T> void PROGRAM<T>::
Change_Instruction(CODE_BLOCK_NODE* N,op_type type,int dest,int src0,int src1)
{
    Remove_Use_Def(N);
    N->inst.type=type;
    N->inst.dest=dest;
    N->inst.src0=src0;
    N->inst.src1=src1;
    Add_Use_Def(N);
}
template struct PROGRAM<float>;
template struct PROGRAM<double>;
}
