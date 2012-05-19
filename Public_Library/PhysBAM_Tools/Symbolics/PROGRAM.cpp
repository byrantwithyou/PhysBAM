#include <PhysBAM_Tools/Data_Structures/HASHTABLE_ITERATOR.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Symbolics/PROGRAM.h>
#include <algorithm>
#include <cstdlib>
#include <cstring>
namespace PhysBAM{

const int mem_mask=0xF0000000;
const int mem_reg=0x00000000;
const int mem_const=0x10000000;
const int mem_in=0x20000000;
const int mem_out=0x30000000;
enum op_flags {flag_none=0,flag_reg_dest=1,flag_reg_src0=2,flag_reg_src1=4};

int op_flags_table[op_last];
int Init_Instructions()
{
    op_flags_table[op_nop]=flag_none;
    op_flags_table[op_copy]=flag_reg_dest|flag_reg_src0;
    op_flags_table[op_add]=flag_reg_dest|flag_reg_src0|flag_reg_src1;
    op_flags_table[op_sub]=flag_reg_dest|flag_reg_src0|flag_reg_src1;
    op_flags_table[op_mul]=flag_reg_dest|flag_reg_src0|flag_reg_src1;
    op_flags_table[op_div]=flag_reg_dest|flag_reg_src0|flag_reg_src1;
    op_flags_table[op_neg]=flag_reg_dest|flag_reg_src0;
    op_flags_table[op_inv]=flag_reg_dest|flag_reg_src0;
    op_flags_table[op_sqrt]=flag_reg_dest|flag_reg_src0;
    op_flags_table[op_exp]=flag_reg_dest|flag_reg_src0;
    op_flags_table[op_ln]=flag_reg_dest|flag_reg_src0;
    op_flags_table[op_pow]=flag_reg_dest|flag_reg_src0|flag_reg_src1;
    op_flags_table[op_lt]=flag_reg_dest|flag_reg_src0|flag_reg_src1;
    op_flags_table[op_le]=flag_reg_dest|flag_reg_src0|flag_reg_src1;
    op_flags_table[op_gt]=flag_reg_dest|flag_reg_src0|flag_reg_src1;
    op_flags_table[op_ge]=flag_reg_dest|flag_reg_src0|flag_reg_src1;
    op_flags_table[op_eq]=flag_reg_dest|flag_reg_src0|flag_reg_src1;
    op_flags_table[op_ne]=flag_reg_dest|flag_reg_src0|flag_reg_src1;
    op_flags_table[op_not]=flag_reg_dest|flag_reg_src0;
    op_flags_table[op_or]=flag_reg_dest|flag_reg_src0|flag_reg_src1;
    op_flags_table[op_and]=flag_reg_dest|flag_reg_src0|flag_reg_src1;
    op_flags_table[op_br_z]=flag_reg_src0;
    op_flags_table[op_br_nz]=flag_reg_src0;
    op_flags_table[op_jmp]=flag_none;
    op_flags_table[op_label]=flag_none;
    return 0;
}
int call_init_instructions=Init_Instructions();
//#####################################################################
// Function Execute
//#####################################################################
template<class T> void PROGRAM<T>::
Execute_Op(ARRAY<T>& reg,int& ip) const
{
    const INSTRUCTION& o=code(ip);
    switch(o.type){
        case op_nop: case op_label: break;
        case op_copy: reg(o.dest)=reg(o.src0);break;
        case op_add: reg(o.dest)=reg(o.src0)+reg(o.src1);break;
        case op_sub: reg(o.dest)=reg(o.src0)-reg(o.src1);break;
        case op_mul: reg(o.dest)=reg(o.src0)*reg(o.src1);break;
        case op_div: reg(o.dest)=reg(o.src0)/reg(o.src1);break;
        case op_neg: reg(o.dest)=-reg(o.src0);break;
        case op_inv: reg(o.dest)=1/reg(o.src0);break;
        case op_sqrt: reg(o.dest)=sqrt(reg(o.src0));break;
        case op_exp: reg(o.dest)=exp(reg(o.src0));break;
        case op_ln: reg(o.dest)=log(reg(o.src0));break;
        case op_pow: reg(o.dest)=pow(reg(o.src0),reg(o.src1));break;
        case op_lt: reg(o.dest)=reg(o.src0)<reg(o.src1);break;
        case op_le: reg(o.dest)=reg(o.src0)<=reg(o.src1);break;
        case op_gt: reg(o.dest)=reg(o.src0)>reg(o.src1);break;
        case op_ge: reg(o.dest)=reg(o.src0)>=reg(o.src1);break;
        case op_eq: reg(o.dest)=reg(o.src0)==reg(o.src1);break;
        case op_ne: reg(o.dest)=reg(o.src0)!=reg(o.src1);break;
        case op_not: reg(o.dest)=!reg(o.src0);break;
        case op_or: reg(o.dest)=reg(o.src0)||reg(o.src1);break;
        case op_and: reg(o.dest)=reg(o.src0)&&reg(o.src1);break;
        case op_br_z: if(!reg(o.src0)) ip=o.dest-1;break;
        case op_br_nz: if(reg(o.src0)) ip=o.dest-1;break;
        case op_jmp: ip=o.dest-1;break;
        default: PHYSBAM_FATAL_ERROR("Missing instruction");
    }
}
//#####################################################################
// Function Execute
//#####################################################################
template<class T> void PROGRAM<T>::
Execute(ARRAY<T>& reg) const
{
    for(int ip=0;ip<code.m;ip++)
        Execute_Op(reg,ip);
}
//#####################################################################
// Function Diff
//#####################################################################
template<class T> int PROGRAM<T>::
Diff(int diff_expr,int diff_var)
{
    ARRAY<INSTRUCTION> c;
    extra_out++;

    ARRAY<int> diff_tmp;
    diff_tmp.Fill(-1);
    for(int i=0;i<var_in.m;i++)
        diff_tmp(dict.Get(var_in(i)))=mem_const;
    diff_tmp(diff_var)=1|mem_const;

    for(int ip=0;ip<code.m;ip++){
        const INSTRUCTION& o=code(ip);
        c.Append(o);
        if(!(op_flags_table[o.type]&flag_reg_dest)) continue;

        int d=-1,s0=-1,s1=-1;
        if(o.dest==diff_expr) continue;
        if(diff_tmp(o.dest)>=0) d=diff_tmp(o.dest);
        else diff_tmp(o.dest)=d=num_tmp++;

        if(op_flags_table[o.type]&flag_reg_src0){
            if((o.src0&mem_mask)==mem_const) s0=mem_const;
            else s0=diff_tmp(o.src0);}

        if(op_flags_table[o.type]&flag_reg_src1){
            if((o.src1&mem_mask)==mem_const) s1=mem_const;
            else s1=diff_tmp(o.src1);}

        switch(o.type){
            case op_copy:case op_neg:{INSTRUCTION in={o.type,d,s0,-1};c.Append(in);}break;
            case op_add:case op_sub:{INSTRUCTION in={o.type,d,s0,s1};c.Append(in);}break;
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
                INSTRUCTION in1={op_mul,d,constants.Append((T).5),in0.dest};c.Append(in1);
                break;}
            case op_exp:{INSTRUCTION in={op_mul,d,s0,o.dest};c.Append(in);}
            case op_ln:{INSTRUCTION in={op_div,d,s0,o.src0};c.Append(in);}
            case op_pow:{
                INSTRUCTION in0={op_ln,num_tmp++,o.src0,-1};c.Append(in0);
                INSTRUCTION in1={op_mul,num_tmp++,s1,in0.dest};c.Append(in1);
                INSTRUCTION in2={op_mul,num_tmp++,o.dest,in1.dest};c.Append(in2);
                INSTRUCTION in3={op_sub,num_tmp++,o.src1,1|mem_const};c.Append(in3);
                INSTRUCTION in4={op_pow,num_tmp++,o.src0,in3.dest};c.Append(in4);
                INSTRUCTION in5={op_mul,num_tmp++,in4.dest,o.src1};c.Append(in5);
                INSTRUCTION in6={op_mul,num_tmp++,in5.dest,s0};c.Append(in6);
                INSTRUCTION in7={op_add,d,in2.dest,in6.dest};c.Append(in7);
                break;}
            case op_lt:case op_le:case op_gt:case op_ge:case op_eq:case op_ne:case op_not:case op_or:case op_and:
                {INSTRUCTION in={op_copy,d,mem_const,-1};c.Append(in);}break;
            default: PHYSBAM_FATAL_ERROR("Missing diff instruction");
        }
    }
    code.Exchange(c);
    return diff_tmp(diff_expr);
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
        binary.Set(">",op_lt);
        binary.Set(">=",op_le);
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
Parse(const char* str)
{
    for(int i=0;i<var_in.m;i++){
        int x=num_tmp++;
        INSTRUCTION in={op_copy,x,i|mem_in,-1};
        code.Append(in);
        dict.Set(var_in(i),x);}
    for(int i=0;i<var_out.m;i++)
        if(!dict.Contains(var_out(i)))
            dict.Set(var_out(i),num_tmp++);

    while(Parse_Command(str)!=-1){}

    for(int i=0;i<var_out.m;i++){
        INSTRUCTION in={op_copy,i|mem_out,dict.Get(var_out(i)),-1};
        code.Append(in);}
}
//#####################################################################
// Function Parse_Command
//#####################################################################
template<class T> int PROGRAM<T>::
Parse_Command(const char*& str)
{
    str+=strspn(str, " \t\n");
    if(!*str || *str==')') return -1;
    puts(str);

    // Number
    char* endptr=0;
    double d=strtod(str,&endptr);
    if(endptr>str){
        str=endptr;
        return constants.Append(d)|mem_const;}

    // Variable
    if(isalpha(*str) || *str=='_'){
        int len=strspn(str, "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_");
        std::string var(str,len);
        str+=len;
        int reg;
        if(dict.Get(var,reg)) return reg;
        reg=dict.Size();
        dict.Set(var,reg);
        return reg;}

    // List
    op_type type;
    if(*str=='('){
        str++;
        str+=strspn(str, " \t\n");
        int len=strcspn(str, " \t\n");
        std::string op(str,len);
        printf("op %s\n",op.c_str());
        str+=len;
        int cur=Parse_Command(str);
        if(op=="if"){
            INSTRUCTION in0={op_br_z,num_labels++,cur,-1};
            code.Append(in0);
            cur=num_tmp++;
            printf("true %s\n", str);
            INSTRUCTION in1={op_copy,cur,Parse_Command(str),-1};
            code.Append(in1);
            INSTRUCTION in2={op_jmp,num_labels++,-1,-1};
            code.Append(in2);
            INSTRUCTION in3={op_label,in0.dest,-1,-1};
            code.Append(in3);
            printf("false %s\n", str);
            INSTRUCTION in4={op_copy,cur,Parse_Command(str),-1};
            code.Append(in4);
            INSTRUCTION in5={op_label,in2.dest,-1,-1};
            code.Append(in5);
            if(Parse_Command(str)!=-1){
                LOG::cout<<"Three arguments expected for 'if'"<<std::endl;
                PHYSBAM_FATAL_ERROR();}}
        else if(cur==-1){
            if(op_tokens.no_arg.Get(op,type)){
                INSTRUCTION in={type,num_tmp++,-1,-1};
                cur=in.dest;
                code.Append(in);}
            else{
                LOG::cout<<"Expected no-arg operator but got '"<<op<<"'"<<std::endl;
                PHYSBAM_FATAL_ERROR();}}
        else{
            int next=Parse_Command(str);
            if(next==-1){
                if(op_tokens.unary.Get(op,type)){
                    INSTRUCTION in={type,num_tmp++,cur,-1};
                    cur=in.dest;
                    code.Append(in);}
                else{
                    LOG::cout<<"Expected unary operator but got '"<<op<<"'"<<std::endl;
                    PHYSBAM_FATAL_ERROR();}}
            else if(op_tokens.chain.Get(op,type)){
                while(next!=-1){
                    INSTRUCTION in={type,num_tmp++,cur,next};
                    cur=in.dest;
                    code.Append(in);
                    next=Parse_Command(str);}}
            else if(op_tokens.binary.Get(op,type)){
                if(Parse_Command(str)!=-1){
                    LOG::cout<<"Expected chaining operator but got '"<<op<<"'"<<std::endl;
                    PHYSBAM_FATAL_ERROR();}
                INSTRUCTION in={type,num_tmp++,cur,next};
                cur=in.dest;
                code.Append(in);}
            else if(op=="=" || op=="setq"){
                if(Parse_Command(str)!=-1){
                    LOG::cout<<"Expected chaining operator but got '"<<op<<"'"<<std::endl;
                    PHYSBAM_FATAL_ERROR();}
                INSTRUCTION in={op_copy,cur,next,-1};
                cur=in.dest;
                code.Append(in);}
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
    "copy  r%d, r%d\n",
    "add   r%d, r%d, r%d\n",
    "sub   r%d, r%d, r%d\n",
    "mul   r%d, r%d, r%d\n",
    "div   r%d, r%d, r%d\n",
    "neg   r%d, r%d\n",
    "inv   r%d, r%d\n",
    "sqrt  r%d, r%d\n",
    "exp   r%d, r%d\n",
    "ln    r%d, r%d\n",
    "pow   r%d, r%d, r%d\n",
    "lt    r%d, r%d, r%d\n",
    "le    r%d, r%d, r%d\n",
    "gt    r%d, r%d, r%d\n",
    "ge    r%d, r%d, r%d\n",
    "eq    r%d, r%d, r%d\n",
    "ne    r%d, r%d, r%d\n",
    "not   r%d, r%d\n",
    "or    r%d, r%d, r%d\n",
    "and   r%d, r%d, r%d\n",
    "br_z  L%d, r%d\n",
    "br_nz L%d, r%d\n",
    "jmp   L%d\n",
    "label L%d\n"
};
//#####################################################################
// Function Print
//#####################################################################
template<class T> void PROGRAM<T>::
Print() const
{
    for(int i=0;i<constants.m;i++)
        printf("const r%d = %g\n", i|mem_const, constants(i));

    for(int i=0;i<code.m;i++){
        printf("% 3d  ",i);
        printf(messages[code(i).type],code(i).dest,code(i).src0,code(i).src1);}
}
//#####################################################################
// Function Finalize
//#####################################################################
template<class T> void PROGRAM<T>::
Finalize()
{
    int off_in=num_tmp,off_out=off_in+var_in.m,off_const=off_out+var_out.m+extra_out;
    ARRAY<int> labels(num_labels);
    int ip2=0;
    for(int ip=0;ip<code.m;ip++){
        INSTRUCTION& o=code(ip);
        if(o.type==op_nop) continue;
        if(o.type==op_label){
            labels(o.dest)=ip2;
            continue;}

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
        code(ip2++)=o;}
    code.Resize(ip2);
    for(int ip=0;ip<code.m;ip++){
        INSTRUCTION& o=code(ip);
        if(o.type==op_br_z || o.type==op_br_nz || o.type==op_jmp)
            o.dest=labels(o.dest);}
}
//#####################################################################
// Function Optimize
//#####################################################################
template<class T> void PROGRAM<T>::
Optimize()
{
    
}
template struct PROGRAM<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template struct PROGRAM<double>;
#endif
}
