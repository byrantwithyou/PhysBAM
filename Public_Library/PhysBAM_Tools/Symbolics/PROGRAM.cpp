#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE_ITERATOR.h>
#include <PhysBAM_Tools/Data_Structures/UNION_FIND.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Symbolics/DOMINATORS.h>
#include <PhysBAM_Tools/Symbolics/PROGRAM.h>
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
        case op_pow: if(in0 && (*in0==0 || *in0==1)){out=*in0;return true;} return false;
            if((in1 && *in1==0)){out=1;return true;} return false;
        case op_or: if((in0 && *in0) || (in1 && *in1)){out=1;return true;} return false;
        case op_and: if((in0 && !*in0) || (in1 && !*in1)){out=0;return true;}} return false;
    return false;
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

    for(int bl=0;bl<code_blocks.m;bl++){
        if(!code_blocks(bl)) continue;
        ARRAY<INSTRUCTION> c;
        ARRAY<INSTRUCTION>& code=code_blocks(bl)->code;
        for(int ip=0;ip<code.m;ip++){
            const INSTRUCTION& o=code(ip);
            c.Append(o);
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
                case op_copy:case op_neg:{INSTRUCTION in={o.type,d,s0,-1};c.Append(in);}break;
                case op_add:case op_sub:case op_phi:{INSTRUCTION in={o.type,d,s0,s1};c.Append(in);}break;
                case op_mul:{
                    INSTRUCTION in0={op_mul,num_tmp++,s0,o.src1};c.Append(in0);
                    INSTRUCTION in1={op_mul,num_tmp++,o.src0,s1};c.Append(in1);
                    INSTRUCTION in2={op_add,d,in0.dest,in1.dest};c.Append(in2);
                    break;}
                case op_div:{
                    INSTRUCTION in0={op_div,num_tmp++,s1,o.src1};c.Append(in0);
                    INSTRUCTION in1={op_mul,num_tmp++,o.dest,in0.dest};c.Append(in1);
                    INSTRUCTION in2={op_div,num_tmp++,s0,o.src1};c.Append(in2);
                    INSTRUCTION in3={op_sub,d,in2.dest,in1.dest};c.Append(in3);
                    break;}
                case op_inv:{
                    INSTRUCTION in0={op_mul,num_tmp++,o.dest,o.dest};c.Append(in0);
                    INSTRUCTION in1={op_mul,num_tmp++,s0,in0.dest};c.Append(in1);
                    INSTRUCTION in2={op_neg,d,in1.dest};c.Append(in2);
                    break;}
                case op_sqrt:{
                    INSTRUCTION in0={op_div,num_tmp++,s0,o.dest};c.Append(in0);
                    INSTRUCTION in1={op_mul,d,Add_Constant((T).5),in0.dest};c.Append(in1);
                    break;}
                case op_exp:{INSTRUCTION in={op_mul,d,s0,o.dest};c.Append(in);}
                case op_ln:{INSTRUCTION in={op_div,d,s0,o.src0};c.Append(in);}
                case op_pow:{
                    INSTRUCTION in0={op_ln,num_tmp++,o.src0,-1};c.Append(in0);
                    INSTRUCTION in1={op_mul,num_tmp++,s1,in0.dest};c.Append(in1);
                    INSTRUCTION in2={op_mul,num_tmp++,o.dest,in1.dest};c.Append(in2);
                    INSTRUCTION in3={op_sub,num_tmp++,o.src1,Add_Constant(1)};c.Append(in3);
                    INSTRUCTION in4={op_pow,num_tmp++,o.src0,in3.dest};c.Append(in4);
                    INSTRUCTION in5={op_mul,num_tmp++,in4.dest,o.src1};c.Append(in5);
                    INSTRUCTION in6={op_mul,num_tmp++,in5.dest,s0};c.Append(in6);
                    INSTRUCTION in7={op_add,d,in2.dest,in6.dest};c.Append(in7);
                    break;}
                case op_lt:case op_le:case op_gt:case op_ge:case op_eq:case op_ne:case op_not:case op_or:case op_and:
                    {INSTRUCTION in={op_copy,d,Add_Constant(0),-1};c.Append(in);}break;
                default: PHYSBAM_FATAL_ERROR("Missing diff instruction");
            }
        }
        code.Exchange(c);
    }
    return var_out.m+extra_out++;
}
struct OPERATOR_DEFINITIONS
{
    HASHTABLE<std::string,op_type> no_arg,unary,binary,chain;

    OPERATOR_DEFINITIONS()
    {
        unary.Set("+",op_copy);
        unary.Set("-",op_neg);
        unary.Set("*",op_copy);
        unary.Set("neg",op_neg);
        unary.Set("inv",op_inv);
        unary.Set("sqrt",op_sqrt);
        unary.Set("exp",op_exp);
        unary.Set("ln",op_ln);
        unary.Set("!",op_not);
        unary.Set("not",op_not);

        binary.Set("+",op_add);
        binary.Set("-",op_sub);
        binary.Set("*",op_mul);
        binary.Set("/",op_div);
        binary.Set("^",op_pow);
        binary.Set("<",op_lt);
        binary.Set("<=",op_le);
        binary.Set(">",op_gt);
        binary.Set(">=",op_ge);
        binary.Set("==",op_eq);
        binary.Set("!=",op_ne);
        binary.Set("||",op_or);
        binary.Set("or",op_or);
        binary.Set("&&",op_and);
        binary.Set("and",op_and);

        chain.Set("+",op_add);
        chain.Set("-",op_sub);
        chain.Set("*",op_mul);
        chain.Set("/",op_div);
        chain.Set("||",op_or);
        chain.Set("or",op_or);
        chain.Set("&&",op_and);
        chain.Set("and",op_and);
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

    CODE_BLOCK<T>* A=new CODE_BLOCK<T>;
    A->id=code_blocks.Append(A);

    for(int i=0;i<var_in.m;i++){
        int x=num_tmp++;
        INSTRUCTION in={op_copy,x,i|mem_in,-1};
        code_blocks(0)->code.Append(in);
        dict.Set(var_in(i),x);}
    for(int i=0;i<var_out.m;i++)
        if(!dict.Contains(var_out(i)))
            dict.Set(var_out(i),num_tmp++);

    LOG::cout<<dict<<"  "<<num_tmp<<std::endl;

    while(Parse_Command(str)!=-1){}

    if(keep_all_vars){
        HASHTABLE<std::string> known;
        known.Set_All(var_out);
        ARRAY<std::string> extra_vars;
        for(HASHTABLE<std::string,int>::ITERATOR it(dict);it.Valid();it.Next())
            if(!known.Contains(it.Key()))
                extra_vars.Append(it.Key());
        extra_vars.Sort();
        var_out.Append_Elements(extra_vars);}

    LOG::cout<<dict<<"  "<<num_tmp<<std::endl;

    LOG::cout<<var_out<<std::endl;
    ARRAY<INSTRUCTION> zero_inst;
    for(int i=0;i<var_out.m;i++){
        INSTRUCTION in1={op_copy,i|mem_out,dict.Get(var_out(i)),-1};
        code_blocks.Last()->code.Append(in1);
        INSTRUCTION in2={op_copy,dict.Get(var_out(i)),Add_Constant(0),-1};
        zero_inst.Append(in2);}
    zero_inst.Append_Elements(code_blocks(0)->code);
    zero_inst.Exchange(code_blocks(0)->code);

    Print();

    Make_SSA();
}
//#####################################################################
// Function Parse_Command
//#####################################################################
template<class T> int PROGRAM<T>::
Parse_Command(const char*& str)
{
    str+=strspn(str, " \t\n");
    if(!*str || *str==')') return -1;

    // Number
    char* endptr=0;
    double d=strtod(str,&endptr);
    if(endptr>str){
        str=endptr;
        return Add_Constant(d);}

    // Variable
    if(isalpha(*str) || *str=='_'){
        int len=strspn(str, "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_");
        std::string var(str,len);
        str+=len;
        int reg;
        if(dict.Get(var,reg)) return reg;
        reg=num_tmp++;
        dict.Set(var,reg);
        return reg;}

    // List
    op_type type;
    if(*str=='('){
        str++;
        str+=strspn(str, " \t\n");
        int len=strcspn(str, " \t\n");
        std::string op(str,len);
        str+=len;
        int cur=Parse_Command(str);
        if(op=="if"){
            CODE_BLOCK<T>* A=code_blocks.Last();
            CODE_BLOCK<T>* B=new CODE_BLOCK<T>;
            CODE_BLOCK<T>* C=new CODE_BLOCK<T>;
            CODE_BLOCK<T>* D=new CODE_BLOCK<T>;
            A->next[0]=B;
            A->next[1]=C;
            D->prev[0]=B;
            D->prev[1]=C;
            B->prev[0]=A;
            B->next[0]=D;
            C->prev[0]=A;
            C->next[0]=D;
            B->id=code_blocks.Append(B);
            int if_out=num_tmp++;
            INSTRUCTION in_B={op_copy,if_out,Parse_Command(str),-1};
            B->code.Append(in_B);
            C->id=code_blocks.Append(C);
            INSTRUCTION in_C={op_copy,if_out,Parse_Command(str),-1};
            C->code.Append(in_C);
            INSTRUCTION in_A={op_br_z,C->id,cur,B->id};
            A->code.Append(in_A);
            cur=if_out;
            D->id=code_blocks.Append(D);
            if(Parse_Command(str)!=-1){
                LOG::cout<<"Three arguments expected for 'if'"<<std::endl;
                PHYSBAM_FATAL_ERROR();}}
        else if(cur==-1){
            if(op_tokens.no_arg.Get(op,type)){
                INSTRUCTION in={type,num_tmp++,-1,-1};
                cur=in.dest;
                code_blocks.Last()->code.Append(in);}
            else{
                LOG::cout<<"Expected no-arg operator but got '"<<op<<"'"<<std::endl;
                PHYSBAM_FATAL_ERROR();}}
        else{
            int next=Parse_Command(str);
            if(next==-1){
                if(op_tokens.unary.Get(op,type)){
                    INSTRUCTION in={type,num_tmp++,cur,-1};
                    cur=in.dest;
                    code_blocks.Last()->code.Append(in);}
                else{
                    LOG::cout<<"Expected unary operator but got '"<<op<<"'"<<std::endl;
                    PHYSBAM_FATAL_ERROR();}}
            else if(op_tokens.chain.Get(op,type)){
                while(next!=-1){
                    INSTRUCTION in={type,num_tmp++,cur,next};
                    cur=in.dest;
                    code_blocks.Last()->code.Append(in);
                    next=Parse_Command(str);}}
            else if(op_tokens.binary.Get(op,type)){
                if(Parse_Command(str)!=-1){
                    LOG::cout<<"Expected chaining operator but got '"<<op<<"'"<<std::endl;
                    PHYSBAM_FATAL_ERROR();}
                INSTRUCTION in={type,num_tmp++,cur,next};
                cur=in.dest;
                code_blocks.Last()->code.Append(in);}
            else if(op=="=" || op=="setq"){
                if(Parse_Command(str)!=-1){
                    LOG::cout<<"Expected chaining operator but got '"<<op<<"'"<<std::endl;
                    PHYSBAM_FATAL_ERROR();}
                INSTRUCTION in={op_copy,cur,next,-1};
                cur=in.dest;
                code_blocks.Last()->code.Append(in);}
            else{
                LOG::cout<<"Unrecognized operator '"<<op<<"'"<<std::endl;
                PHYSBAM_FATAL_ERROR();}}
        str+=strspn(str, " \t\n");
        if(*str++==')') return cur;
        LOG::cout<<"Expected ) but got '"<<*--str<<"'"<<std::endl;
        PHYSBAM_FATAL_ERROR();}

    LOG::cout<<"Expected token but got '"<<*str<<"'"<<std::endl;
    PHYSBAM_FATAL_ERROR();
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
Print(const ARRAY<INSTRUCTION>& code) const
{
    for(int i=0;i<code.m;i++){
        printf("% 3d  ",i);
        printf(messages[code(i).type],"rcio"[code(i).dest>>mem_shift],code(i).dest&~mem_mask,
            "rcio"[code(i).src0>>mem_shift],code(i).src0&~mem_mask,
            "rcio"[code(i).src1>>mem_shift],code(i).src1&~mem_mask);}
}
//#####################################################################
// Function Print
//#####################################################################
template<class T> void PROGRAM<T>::
Print() const
{
    for(int i=0;i<constants.m;i++)
        printf("const c%d = %g\n", i, constants(i));

    printf("block code:\n");
    for(int i=0;i<code_blocks.m;i++){
        if(!code_blocks(i)) continue;
        if(code_blocks(i)->prev[0])
            printf(messages[op_label],i);
        Print(code_blocks(i)->code);
        if(code_blocks(i)->next[0])
            printf(messages[op_jmp],code_blocks(i)->next[0]->id);}

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
    Leave_SSA();

    int off_in=num_tmp,off_out=off_in+var_in.m,off_const=off_out+var_out.m+extra_out;
    ARRAY<int> labels(code_blocks.m);
    flat_code.Remove_All();
    for(int bl=0;bl<code_blocks.m;bl++){
        if(!code_blocks(bl)) continue;
        labels(bl)=flat_code.m;
        if(flat_code.m && flat_code.Last().type==op_jmp && flat_code.Last().dest==bl)
            flat_code.Pop();
        ARRAY<INSTRUCTION>& code=code_blocks(bl)->code;
        for(int ip=0;ip<code.m;ip++){
            INSTRUCTION o=code(ip);
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
        if(code.m && code_blocks(bl)->next[0] && code.Last().type!=op_br_z && code.Last().type!=op_br_nz && code.Last().type!=op_jmp){
            INSTRUCTION in={op_jmp,code_blocks(bl)->next[0]->id,0,0};
            flat_code.Append(in);}}
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
    Local_Common_Expresssion_Elimination();
    Remove_Dead_Code();
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

    LOG::cout<<dom.idom<<std::endl;

    LOG::cout<<graph_out<<std::endl;

    LOG::cout<<dom.frontier<<std::endl;

    ARRAY<int> has_already(code_blocks.m),work(code_blocks.m),V;
    has_already.Fill(-1);
    work.Fill(-1);
    HASHTABLE<int,HASHTABLE<int> > HA;
    V=IDENTITY_ARRAY<>(num_tmp);
    HASHTABLE<int> V_used;
    V_used.Set_All(V);
    for(int bl=0;bl<code_blocks.m;bl++){
        if(!code_blocks(bl)) continue;
        ARRAY<INSTRUCTION>& code=code_blocks(bl)->code;
        for(int ip=0;ip<code.m;ip++){
            const INSTRUCTION& o=code(ip);
            if(op_flags_table[o.type]&flag_reg_dest)
                HA.Get_Or_Insert(o.dest).Set(bl);
            if(op_flags_table[o.type]&flag_reg_src0)
                HA.Get_Or_Insert(o.src0).Set(bl);
            if(op_flags_table[o.type]&flag_reg_src1)
                HA.Get_Or_Insert(o.src1).Set(bl);}}

    ARRAY<ARRAY<INSTRUCTION> > new_phi(code_blocks.m);

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
                        new_phi(y).Append(in);}
                    has_already(y)=vi;
                    if(work(y)<vi){
                        work(y)=vi;
                        W.Append(y);}}}}}

    for(int i=0;i<new_phi.m;i++)
        if(new_phi(i).m){
            new_phi(i).Append_Elements(code_blocks(i)->code);
            new_phi(i).Exchange(code_blocks(i)->code);}

    ARRAY<ARRAY<int> > S(V.m);
    int undef=Add_Constant(0);
    for(int i=0;i<S.m;i++)
        S(i).Append(undef);
    Make_SSA_Relabel(num_tmp,S,V_used,dom.children,0);
    for(int i=0;i<S.m;i++){
        assert(S(i).m==1 && S(i)(0)==undef);}

    LOG::cout<<__FUNCTION__<<std::endl;
    Print();
}
//#####################################################################
// Function Make_SSA_Relabel
//#####################################################################
template<class T> void PROGRAM<T>::
Make_SSA_Relabel(int& count,ARRAY<ARRAY<int> >& S,HASHTABLE<int>& V_used,const ARRAY<ARRAY<int> >& children,int bl)
{
    LOG::cout<<count<<std::endl;
    ARRAY<int> S_pop;
    for(int ip=0;ip<code_blocks(bl)->code.m;ip++){
        INSTRUCTION& o=code_blocks(bl)->code(ip);
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
        if(CODE_BLOCK<T>* nb=code_blocks(bl)->next[s]){
            int j=nb->prev[1]==code_blocks(bl);
            for(int i=0;i<nb->code.m && nb->code(i).type==op_phi;i++){
                int* reg=j?&nb->code(i).src1:&nb->code(i).src0;
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
    for(int i=0;i<code_blocks.m;i++){
        if(!code_blocks(i)) continue;
        for(int j=0;j<code_blocks(i)->code.m;j++){
            INSTRUCTION& o=code_blocks(i)->code(j);
            if(o.type!=op_phi) break;
            if(o.src0==o.src1) var_map(o.dest)=o.src0;
            else{
                INSTRUCTION in0={op_copy,o.dest,o.src0,-1};
                INSTRUCTION in1={op_copy,o.dest,o.src1,-1};
                code_blocks(i)->prev[0]->code.Append(in0);
                code_blocks(i)->prev[1]->code.Append(in1);}
            o.type=op_nop;}}
    Relabel_Registers(var_map);

    LOG::cout<<__FUNCTION__<<std::endl;
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
    defs(0)=CONSTANT_ARRAY<VECTOR<int,2> >(num_tmp,VECTOR<int,2>(-1,-1));
    defs(3)=CONSTANT_ARRAY<VECTOR<int,2> >(var_out.m+extra_out,VECTOR<int,2>(-1,-1));

    for(int bl=0;bl<code_blocks.m;bl++){
        if(!code_blocks(bl)) continue;
        for(int ip=0;ip<code_blocks(bl)->code.m;ip++){
            INSTRUCTION& o=code_blocks(bl)->code(ip);
            VECTOR<int,2> v(bl,ip);
            if(op_flags_table[o.type]&flag_reg_dest) Set_Def(o.dest,v);
            if(op_flags_table[o.type]&flag_reg_src0) Add_Use(o.src0,v);
            if(op_flags_table[o.type]&flag_reg_src1 && o.src0!=o.src1) Add_Use(o.src1,v);}}
    def_use_stale=false;
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
    LOG::cout<<new_map<<std::endl;

    for(int bl=0;bl<code_blocks.m;bl++){
        if(!code_blocks(bl)) continue;
        for(int ip=0;ip<code_blocks(bl)->code.m;ip++){
            INSTRUCTION& o=code_blocks(bl)->code(ip);
            if(op_flags_table[o.type]&flag_reg_dest && (o.dest&mem_mask)==mem_reg)
                o.dest=new_map(o.dest);
            if(op_flags_table[o.type]&flag_reg_src0 && (o.src0&mem_mask)==mem_reg)
                o.src0=new_map(o.src0);
            if(op_flags_table[o.type]&flag_reg_src1 && (o.src1&mem_mask)==mem_reg)
                o.src1=new_map(o.src1);}}

    num_tmp=k;
    for(HASHTABLE<std::string,int>::ITERATOR it(dict);it.Valid();it.Next())
        it.Data()=new_map(it.Data());

    var_map.Exchange(new_map);
}
//#####################################################################
// Function SCCP_Visit_Instruction
//#####################################################################
template<class T> void PROGRAM<T>::
SCCP_Visit_Instruction(const VECTOR<int,2>& I,VECTOR<ARRAY<VARIABLE_STATE>,4>& variable_state,
    ARRAY<CODE_BLOCK<T>*>& block_worklist,ARRAY<VECTOR<int,2> >& op_worklist,ARRAY<bool>& block_exec)
{
    INSTRUCTION& o=code_blocks(I.x)->code(I.y);
    if(o.type==op_phi){
        VARIABLE_STATE& d=variable_state(o.dest>>mem_shift)(o.dest&~mem_mask);
        if(d.state==VARIABLE_STATE::overdetermined) return;
        VARIABLE_STATE s0=variable_state(o.src0>>mem_shift)(o.src0&~mem_mask);
        VARIABLE_STATE s1=variable_state(o.src1>>mem_shift)(o.src1&~mem_mask);

        if(s0.state==VARIABLE_STATE::overdetermined || s1.state==VARIABLE_STATE::overdetermined){
            d.state=VARIABLE_STATE::overdetermined;
            op_worklist.Append(I);
            return;}

        if(s0.state==VARIABLE_STATE::constant && s1.state==VARIABLE_STATE::constant){
            if(s0.value==s1.value) d.state=VARIABLE_STATE::overdetermined;
            else if(d.state==VARIABLE_STATE::constant) return;
            else{d.state=VARIABLE_STATE::constant;d.value=s0.value;}
            op_worklist.Append(I);
            return;}

        if(s0.state==VARIABLE_STATE::constant || s1.state==VARIABLE_STATE::constant){
            if(d.state==VARIABLE_STATE::constant) return;
            d.state=VARIABLE_STATE::constant;
            d.value=s0.state==VARIABLE_STATE::constant?s0.value:s1.value;
            op_worklist.Append(I);
            return;}
        return;}

    if(op_flags_table[o.type]&flag_reg_dest){
        VARIABLE_STATE& d=variable_state(o.dest>>mem_shift)(o.dest&~mem_mask);
        if(d.state==VARIABLE_STATE::constant || d.state==VARIABLE_STATE::overdetermined) return;
        VARIABLE_STATE* s0=op_flags_table[o.type]&flag_reg_src0?&variable_state(o.src0>>mem_shift)(o.src0&~mem_mask):0;
        VARIABLE_STATE* s1=op_flags_table[o.type]&flag_reg_src1?&variable_state(o.src1>>mem_shift)(o.src1&~mem_mask):0;

        // rhs is known; evaluate instruction
        if((!s0 || s0->state==VARIABLE_STATE::constant) && (!s1 || s1->state==VARIABLE_STATE::constant)){
            d.state=VARIABLE_STATE::constant;
            d.value=Evaluate_Op(o.type,s0?s0->value:0,s1?s1->value:0);
            op_worklist.Append(I);
            return;}

        // Might be able to deduce result with partial information
        if(op_flags_table[o.type]&flag_can_deduce){
            // rhs is partially known; see if result can be deduced anyway
            if((s0 && s0->state==VARIABLE_STATE::constant) || (s1 && s1->state==VARIABLE_STATE::constant)){
                T* a=(s0 && s0->state==VARIABLE_STATE::constant)?&s0->value:0;
                T* b=(s1 && s1->state==VARIABLE_STATE::constant)?&s1->value:0;
                if(Deduce_Op(o.type,d.value,a,b)){
                    d.state=VARIABLE_STATE::constant;
                    op_worklist.Append(I);
                    return;}}
            if((!s0 || s0->state==VARIABLE_STATE::overdetermined) && (!s1 || s1->state==VARIABLE_STATE::overdetermined)){
                d.state=VARIABLE_STATE::overdetermined;
                op_worklist.Append(I);
                return;}}
        else if((s0 && s0->state==VARIABLE_STATE::overdetermined) || (s1 && s1->state==VARIABLE_STATE::overdetermined)){
            d.state=VARIABLE_STATE::overdetermined;
            op_worklist.Append(I);
            return;}
        return;}

    if(o.type==op_nop || o.type==op_label || o.type==op_jmp) return;
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

    if((o.type==op_br_nz)==(s0.value!=0)){
        o.dest=o.src1;
        exchange(code_blocks(I.x)->next[0],code_blocks(I.x)->next[1]);}
    if(!block_exec(o.dest)){
        block_exec(o.dest)=true;
        block_worklist.Append(code_blocks(o.dest));}

    o.type=op_nop;
    Detach_Next(code_blocks(I.x),1);
}
//#####################################################################
// Function Sparse_Conditional_Constant_Propagation
//#####################################################################
template<class T> void PROGRAM<T>::
Sparse_Conditional_Constant_Propagation()
{
    if(def_use_stale) Update_Use_Def();

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

    ARRAY<CODE_BLOCK<T>*> block_worklist;
    ARRAY<VECTOR<int,2> > op_worklist;
    ARRAY<bool> block_exec(code_blocks.m);
    block_exec(0)=true;
    block_worklist.Append(code_blocks(0));

    while(block_worklist.m || op_worklist.m){
        while(op_worklist.m){
            VECTOR<int,2> I=op_worklist.Pop();
            for(HASHTABLE<VECTOR<int,2> >::CONST_ITERATOR it(Get_Uses(I));it.Valid();it.Next()){
                VECTOR<int,2> J=it.Key();
                if(block_exec(J.x))
                    SCCP_Visit_Instruction(J,variable_state,block_worklist,op_worklist,block_exec);}}

        while(block_worklist.m){
            CODE_BLOCK<T>* block=block_worklist.Pop();
            for(int ip=0;ip<block->code.m;ip++)
                SCCP_Visit_Instruction(VECTOR<int,2>(block->id,ip),variable_state,block_worklist,op_worklist,block_exec);}}

    for(int i=0;i<4;i++)
        for(int j=0;j<defs(i).m;j++){
            if(variable_state(i)(j).state==VARIABLE_STATE::constant){
                VECTOR<int,2> I=defs(i)(j);
                if(I.x==-1) continue;
                INSTRUCTION& o=Get_Instruction(I);
                o.type=op_copy;
                LOG::cout<<"solved: "<<I<<"   r"<<j<<" = "<<variable_state(i)(j).value<<std::endl;
                o.src0=Add_Constant(variable_state(i)(j).value);}}
    def_use_stale=true;

    LOG::cout<<__FUNCTION__<<std::endl;
    Print();
}
//#####################################################################
// Function Copy_Propagation
//#####################################################################
template<class T> void PROGRAM<T>::
Propagate_Copy(int old_var,int new_var)
{
    HASHTABLE<VECTOR<int,2> >& old_uses = uses(old_var>>mem_shift)(old_var&~mem_mask);
    HASHTABLE<VECTOR<int,2> >& new_uses = uses(new_var>>mem_shift)(new_var&~mem_mask);
    for(HASHTABLE<VECTOR<int,2> >::ITERATOR it(old_uses);it.Valid();it.Next()){
        VECTOR<int,2> I=it.Key();
        INSTRUCTION& o2=code_blocks(I.x)->code(I.y);
        if(o2.src0==old_var) o2.src0=new_var;
        if(o2.src1==old_var) o2.src1=new_var;
        new_uses.Set(I);}
    old_uses.Remove_All();
}
//#####################################################################
// Function Copy_Propagation
//#####################################################################
template<class T> void PROGRAM<T>::
Copy_Propagation()
{
    if(def_use_stale) Update_Use_Def();

    // Copy propagation
    for(int bl=0;bl<code_blocks.m;bl++){
        if(!code_blocks(bl)) continue;
        for(int ip=0;ip<code_blocks(bl)->code.m;ip++){
            INSTRUCTION& o=code_blocks(bl)->code(ip);
            if(o.type==op_copy)
                Propagate_Copy(o.dest,o.src0);}}

    LOG::cout<<__FUNCTION__<<std::endl;
    Print();
}
//#####################################################################
// Function Simplify_Phis
//#####################################################################
template<class T> void PROGRAM<T>::
Simplify_Phis()
{
    // Replace phi's with unreachble input by a copy.
    for(int bl=0;bl<code_blocks.m;bl++){
        if(!code_blocks(bl)) continue;
        if(!code_blocks(bl)->prev[0] || !code_blocks(bl)->prev[1]){
            int j=!code_blocks(bl)->prev[0];
            for(int ip=0;ip<code_blocks(bl)->code.m;ip++){
                INSTRUCTION& o=code_blocks(bl)->code(ip);
                if(o.type!=op_phi) break;
                o.type=op_copy;
                if(j) o.src0=o.src1;}
            code_blocks(bl)->prev[0]=code_blocks(bl)->prev[j];}}
    def_use_stale=true;

    LOG::cout<<__FUNCTION__<<std::endl;
    Print();
}
//#####################################################################
// Function Remove_Dead_Code
//#####################################################################
template<class T> void PROGRAM<T>::
Remove_Dead_Code()
{
        // Prune unreachable basic blocks
    for(int bl=0;bl<code_blocks.m;bl++){
        CODE_BLOCK<T>* B=code_blocks(bl);
        if(!B) continue;
        if(bl>0 && !B->prev[0] && !B->prev[1]){
            for(int s=0;s<2;s++) if(B->next[s]) Detach_Next(B,s);
            delete B;
            code_blocks(bl)=0;}}

    // Merge basic blocks
    for(int bl=0;bl<code_blocks.m;bl++){
        CODE_BLOCK<T>* B=code_blocks(bl);
        if(!B) continue;
        while(B->next[0] && !B->next[1] && !B->next[0]->prev[1]){
            CODE_BLOCK<T>* nb=B->next[0];
            for(int s=0;s<2;s++){
                if(nb->next[s]){
                    int j=nb->next[s]->prev[1]==nb;
                    nb->next[s]->prev[j]=B;}
                B->next[s]=nb->next[s];}
            B->code.Append_Elements(nb->code);
            delete code_blocks(nb->id);
            code_blocks(nb->id)=0;}}

    Update_Use_Def();
    ARRAY<int> worklist;
    for(int j=0;j<defs(0).m;j++)
        if(defs(0)(j).x!=-1 && !uses(0)(j).Size())
            worklist.Append(j);

    while(worklist.m){
        int var=worklist.Pop();
        VECTOR<int,2> I=Get_Def(var);
        INSTRUCTION& o=Get_Instruction(I);
        if(op_flags_table[o.type]&flag_reg_src0 && (o.src0&mem_mask)==mem_reg){
            HASHTABLE<VECTOR<int,2> >& h=Get_Uses(o.src0);
            if(h.Delete_If_Present(I) && !h.Size())
                worklist.Append(o.src0);}
        if(op_flags_table[o.type]&flag_reg_src1 && (o.src1&mem_mask)==mem_reg){
            HASHTABLE<VECTOR<int,2> >& h=Get_Uses(o.src1);
            if(h.Delete_If_Present(I) && !h.Size())
                worklist.Append(o.src1);}
        o.type=op_nop;
        Set_Def(var,VECTOR<int,2>(-1,-1));}

    for(int bl=0;bl<code_blocks.m;bl++){
        if(!code_blocks(bl)) continue;
        int k=0;
        for(int ip=0;ip<code_blocks(bl)->code.m;ip++){
            INSTRUCTION o=code_blocks(bl)->code(ip);
            if(o.type!=op_label && o.type!=op_nop)
                code_blocks(bl)->code(k++)=o;}
        code_blocks(bl)->code.Resize(k);}

    def_use_stale=true;

    LOG::cout<<__FUNCTION__<<std::endl;
    Print();
}
//#####################################################################
// Function Detach_Next
//#####################################################################
template<class T> void PROGRAM<T>::
Detach_Next(CODE_BLOCK<T>* B,int j)
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
    constant_lookup[x]=i;
    return i|mem_const;
}
//#####################################################################
// Function Local_Common_Expresssion_Elimination
//#####################################################################
template<class T> void PROGRAM<T>::
Local_Common_Expresssion_Elimination(CODE_BLOCK<T>* B)
{
    ARRAY<INSTRUCTION>& code=B->code;
    HASHTABLE<VECTOR<int,3>,int> expr;
    for(int ip=0;ip<code.m;ip++){
        INSTRUCTION& o=code(ip);
        if(!(op_flags_table[o.type]&flag_reg_dest)) continue;
        Reduce_In_Place(o);
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
    for(int bl=0;bl<code_blocks.m;bl++){
        if(!code_blocks(bl)) continue;
        Local_Common_Expresssion_Elimination(code_blocks(bl));}
}
//#####################################################################
// Function Reduce_In_Place
//#####################################################################
template<class T> void PROGRAM<T>::
Reduce_In_Place(INSTRUCTION& o)
{
    int const_0=Add_Constant(0),const_1=Add_Constant(1);
    switch(o.type){
        case op_add:
            if(o.src1==const_0){o.type=op_copy;return;}
            if(o.src0==const_0){o.type=op_copy;o.src0=o.src1;return;}
            break;
        case op_sub:
            if(o.src0==o.src1){o.type=op_copy;o.src0=const_0;return;}
            if(o.src1==const_0){o.type=op_copy;return;}
            if(o.src0==const_0){o.type=op_neg;o.src0=o.src1;return;}
            break;
        case op_mul:
            if(o.src1==const_0 || o.src0==const_0){o.type=op_copy;o.src0=const_0;return;}
            if(o.src1==const_1){o.type=op_copy;return;}
            if(o.src0==const_1){o.type=op_copy;o.src0=o.src1;return;}
            break;
        case op_div:
            if(o.src0==o.src1){o.type=op_copy;o.src0=const_1;return;}
            if(o.src0==const_0 || o.src1==const_1){o.type=op_copy;return;}
            if(o.src0==const_1){o.type=op_inv;o.src0=o.src1;return;}
            break;
        case op_neg: break;
        case op_inv: break;
        case op_sqrt: break;
        case op_exp: break;
        case op_ln: break;
        case op_pow:
            if(o.src1==const_0){o.type=op_copy;o.src0=const_0;return;}
            if(o.src1==const_1){o.type=op_copy;return;}
            if(o.src1==Add_Constant(2)){o.type=op_mul;o.src1=o.src0;return;}
            if(o.src1==Add_Constant(-1)){o.type=op_neg;return;}
            break;
        case op_lt:
        case op_gt:
        case op_ne:
            if(o.src0==o.src1){o.type=op_copy;o.src0=const_0;return;}
            break;
        case op_le:
        case op_ge:
        case op_eq:
            if(o.src0==o.src1){o.type=op_copy;o.src0=const_1;return;}
            break;
        case op_or:
        case op_and:
            if(o.src0==o.src1){o.type=op_copy;return;}
            break;
        default: break;}
}
template struct PROGRAM<float>;
template struct PROGRAM<double>;
}
