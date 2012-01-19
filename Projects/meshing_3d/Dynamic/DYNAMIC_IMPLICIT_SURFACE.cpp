#include <PhysBAM_Tools/Data_Structures/HASHTABLE_ITERATOR.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <algorithm>
#include <cstdlib>
#include <cstring>
#include "DYNAMIC_IMPLICIT_SURFACE.h"
namespace PhysBAM{

template<class T>
struct DYNAMIC_COORD:public DYNAMIC_OP<T>
{
    typedef VECTOR<T,3> TV;using DYNAMIC_OP<T>::S;using DYNAMIC_OP<T>::D;using DYNAMIC_OP<T>::H;using DYNAMIC_OP<T>::A; \
    int c;
    DYNAMIC_COORD(): c(0) {}
    virtual DYNAMIC_OP<T>* New() const {return new DYNAMIC_COORD<T>;}
    virtual ~DYNAMIC_COORD() {}
    virtual void Compute(const TV& X) {PHYSBAM_ASSERT(A.m==0 && c>=1 && c<=3);S=X(c);D=TV::Axis_Vector(c);}
};

#define DYNAMIC_TEMPLATE(NAME) \
template<class T> \
struct NAME:public DYNAMIC_OP<T> \
{ \
    typedef VECTOR<T,3> TV;using DYNAMIC_OP<T>::S;using DYNAMIC_OP<T>::D;using DYNAMIC_OP<T>::H;using DYNAMIC_OP<T>::A; \
    virtual ~NAME() {} \
    virtual DYNAMIC_OP<T>* New() const {return new NAME;} \
    virtual void Compute(const TV& X); \
}

DYNAMIC_TEMPLATE(DYNAMIC_NUMBER);
template<class T> void DYNAMIC_NUMBER<T>::Compute(const TV& X)
{
    PHYSBAM_ASSERT(A.m==0);
}

DYNAMIC_TEMPLATE(DYNAMIC_IF);
template<class T> void DYNAMIC_IF<T>::Compute(const TV& X)
{
    PHYSBAM_ASSERT(A.m==3);
    DYNAMIC_OP<T> *P=A(1)->S?A(2):A(3);
    S=P->S;
    D=P->D;
    H=P->H;
}

DYNAMIC_TEMPLATE(DYNAMIC_ABS);
template<class T> void DYNAMIC_ABS<T>::Compute(const TV& X)
{
    PHYSBAM_ASSERT(A.m==1);
    if(A(1)->S>=0){S=A(1)->S;D=A(1)->D;H=A(1)->H;}
    else{S=-A(1)->S;D=-A(1)->D;H=-A(1)->H;}
}

DYNAMIC_TEMPLATE(DYNAMIC_ADD);
template<class T> void DYNAMIC_ADD<T>::Compute(const TV& X)
{
    PHYSBAM_ASSERT(A.m>=1);
    S=A(1)->S;
    D=A(1)->D;
    H=A(1)->H;
    for(int j=2;j<=A.m;j++){S+=A(j)->S;D+=A(j)->D;H+=A(j)->H;}
}

DYNAMIC_TEMPLATE(DYNAMIC_SUB);
template<class T> void DYNAMIC_SUB<T>::Compute(const TV& X)
{
    PHYSBAM_ASSERT(A.m>=1);
    if(A.m==1){
        S=-A(1)->S;
        D=-A(1)->D;
        H=-A(1)->H;}
    else{
        S=A(1)->S;
        D=A(1)->D;
        H=A(1)->H;
        for(int j=2;j<=A.m;j++){S-=A(j)->S;D-=A(j)->D;H-=A(j)->H;}}
}

DYNAMIC_TEMPLATE(DYNAMIC_MUL);
template<class T> void DYNAMIC_MUL<T>::Compute(const TV& X)
{
    PHYSBAM_ASSERT(A.m>=1);
    S=A(1)->S;
    D=A(1)->D;
    H=A(1)->H;
    for(int j=2;j<=A.m;j++){
        H=A(j)->S*H+S*A(j)->H+MATRIX<T,3>::Outer_Product(A(j)->D,D)+MATRIX<T,3>::Outer_Product(D,A(j)->D);
        D=S*A(j)->D+A(j)->S*D;
        S*=A(j)->S;}
}

DYNAMIC_TEMPLATE(DYNAMIC_DIV);
template<class T> void DYNAMIC_DIV<T>::Compute(const TV& X)
{
    PHYSBAM_ASSERT(A.m>=1);
    S=A(1)->S;
    D=A(1)->D;
    H=A(1)->H;
    for(int j=2;j<=A.m;j++){
        T SI=1/A(j)->S;
        TV DI=-SI*SI*A(j)->D;
        MATRIX<T,3> HI=-SI*SI*A(j)->H-2*SI*MATRIX<T,3>::Outer_Product(DI,A(j)->D);
        H=SI*H+S*HI+MATRIX<T,3>::Outer_Product(DI,D)+MATRIX<T,3>::Outer_Product(D,DI);
        D=S*DI+SI*D;
        S*=SI;}
}

DYNAMIC_TEMPLATE(DYNAMIC_RECIP);
template<class T> void DYNAMIC_RECIP<T>::Compute(const TV& X)
{
    PHYSBAM_ASSERT(A.m==1);
    S=1/A(1)->S;
    D=-S*S*A(1)->D;
    H=-S*S*A(1)->H-2*S*MATRIX<T,3>::Outer_Product(D,A(1)->D);
}

DYNAMIC_TEMPLATE(DYNAMIC_SQRT);
template<class T> void DYNAMIC_SQRT<T>::Compute(const TV& X)
{
    PHYSBAM_ASSERT(A.m==1);
    if(A(1)->S<=0){S=0;D=TV();H=MATRIX<T,3>();}
    else{
        S=sqrt(A(1)->S);
        D=A(1)->D/(2*S);
        H=A(1)->H/(2*S)-MATRIX<T,3>::Outer_Product(D,A(1)->D);}
}

DYNAMIC_TEMPLATE(DYNAMIC_LT);
template<class T> void DYNAMIC_LT<T>::Compute(const TV& X)
{
    PHYSBAM_ASSERT(A.m==2);
    S=A(1)->S<A(2)->S;
}

DYNAMIC_TEMPLATE(DYNAMIC_GT);
template<class T> void DYNAMIC_GT<T>::Compute(const TV& X)
{
    PHYSBAM_ASSERT(A.m==2);
    S=A(1)->S>A(2)->S;
}

DYNAMIC_TEMPLATE(DYNAMIC_LE);
template<class T> void DYNAMIC_LE<T>::Compute(const TV& X)
{
    PHYSBAM_ASSERT(A.m==2);
    S=A(1)->S<=A(2)->S;
}

DYNAMIC_TEMPLATE(DYNAMIC_GE);
template<class T> void DYNAMIC_GE<T>::Compute(const TV& X)
{
    PHYSBAM_ASSERT(A.m==2);
    S=A(1)->S>=A(2)->S;
}

DYNAMIC_TEMPLATE(DYNAMIC_NE);
template<class T> void DYNAMIC_NE<T>::Compute(const TV& X)
{
    PHYSBAM_ASSERT(A.m==2);
    S=A(1)->S!=A(2)->S;
}

DYNAMIC_TEMPLATE(DYNAMIC_EQ);
template<class T> void DYNAMIC_EQ<T>::Compute(const TV& X)
{
    PHYSBAM_ASSERT(A.m==2);
    S=A(1)->S==A(2)->S;
}

DYNAMIC_TEMPLATE(DYNAMIC_NOT);
template<class T> void DYNAMIC_NOT<T>::Compute(const TV& X)
{
    PHYSBAM_ASSERT(A.m==2);
    S=!A(1)->S;
}

DYNAMIC_TEMPLATE(DYNAMIC_AND);
template<class T> void DYNAMIC_AND<T>::Compute(const TV& X)
{
    PHYSBAM_ASSERT(A.m==2);
    S=(A(1)->S&&A(2)->S);
}

DYNAMIC_TEMPLATE(DYNAMIC_OR);
template<class T> void DYNAMIC_OR<T>::Compute(const TV& X)
{
    PHYSBAM_ASSERT(A.m==2);
    S=(A(1)->S||A(2)->S);
}
//#####################################################################
// Destructor
//#####################################################################
template<class T> DYNAMIC_IMPLICIT_SURFACE<T>::
~DYNAMIC_IMPLICIT_SURFACE()
{
    eval_list.Delete_Pointers_And_Clean_Memory();
    dict.Delete_Pointers_Stored_In_Table();
    token_lookup.Delete_Pointers_Stored_In_Table();
}
//#####################################################################
// Function Initialize_Tokens
//#####################################################################
template<class T> void DYNAMIC_IMPLICIT_SURFACE<T>::
Initialize_Tokens()
{
    token_lookup.Set("if",new DYNAMIC_IF<T>);
    token_lookup.Set("+",new DYNAMIC_ADD<T>);
    token_lookup.Set("-",new DYNAMIC_SUB<T>);
    token_lookup.Set("*",new DYNAMIC_MUL<T>);
    token_lookup.Set("/",new DYNAMIC_DIV<T>);
    token_lookup.Set("inv",new DYNAMIC_RECIP<T>);
    token_lookup.Set("<",new DYNAMIC_LT<T>);
    token_lookup.Set(">",new DYNAMIC_GT<T>);
    token_lookup.Set("<=",new DYNAMIC_LE<T>);
    token_lookup.Set(">=",new DYNAMIC_GE<T>);
    token_lookup.Set("!=",new DYNAMIC_NE<T>);
    token_lookup.Set("==",new DYNAMIC_EQ<T>);
    token_lookup.Set("!",new DYNAMIC_NOT<T>);
    token_lookup.Set("or",new DYNAMIC_OR<T>);
    token_lookup.Set("and",new DYNAMIC_AND<T>);
    token_lookup.Set("abs",new DYNAMIC_ABS<T>);
    token_lookup.Set("sqrt",new DYNAMIC_SQRT<T>);
}
//#####################################################################
// Function Next_Token
//#####################################################################
template<class T> std::string DYNAMIC_IMPLICIT_SURFACE<T>::
Next_Token(const char*& str) const
{
    str+=strspn(str, " \t\n");
    if(*str=='('){str++;return "(";}
    if(*str==')'){str++;return ")";}
    int len=strcspn(str, " \t\n()");
    std::string s(str,len);
    str+=len;
    return s;
}
//#####################################################################
// Function Parse
//#####################################################################
template<class T> DYNAMIC_OP<T>* DYNAMIC_IMPLICIT_SURFACE<T>::
Parse(const char*& str)
{
    const char* orig=str;
    std::string token=Next_Token(str);
    if(token.empty()) return 0;
    if(token==")"){str=orig;return 0;}

    // Number
    char* endptr=0;
    double d=strtod(token.c_str(),&endptr);
    if(!*endptr){
        DYNAMIC_NUMBER<T>* op=new DYNAMIC_NUMBER<T>;
        op->S=d;
        eval_list.Append(op);
        return op;}

    // Variable
    if(DYNAMIC_OP<T>** var=dict.Get_Pointer(token)) return *var;

    // List
    if(token=="("){
        std::string next=Next_Token(str);
        if(next=="setq"){
            std::string name=Next_Token(str);
            DYNAMIC_OP<T>* val=Parse(str);
            std::string paren=Next_Token(str);
            if(paren!=")"){
                LOG::cout<<"Found token '"<<paren<<"' but expected ')'"<<std::endl;
                PHYSBAM_FATAL_ERROR();}
            dict.Set(name,val);
            return val;}

        DYNAMIC_OP<T>** pop=token_lookup.Get_Pointer(next);
        if(!pop){
            LOG::cout<<"Found token '"<<next<<"' but expected operator"<<std::endl;
            PHYSBAM_FATAL_ERROR();}
        DYNAMIC_OP<T>* op=(*pop)->New();

        while(DYNAMIC_OP<T>* arg=Parse(str)) op->A.Append(arg);

        next=Next_Token(str);
        if(next!=")"){
            LOG::cout<<"Found token '"<<next<<"' but expected ')'"<<std::endl;
            PHYSBAM_FATAL_ERROR();}
        eval_list.Append(op);
        return op;}

    // Coordinate
    if(token=="x" || token=="y" || token=="z"){
        DYNAMIC_COORD<T>* op=new DYNAMIC_COORD<T>;
        op->c=token[0]-'x'+1;
        eval_list.Append(op);
        return op;}

    LOG::cout<<"Unexpected token '"<<token<<std::endl;
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function Eval
//#####################################################################
template<class T> void DYNAMIC_IMPLICIT_SURFACE<T>::
Eval(const TV& X) const
{
    for(int i=0;i<eval_list.m;i++) eval_list(i)->Compute(X);

//    for(typename HASHTABLE<std::string,DYNAMIC_OP<T>*>::CONST_ITERATOR it(dict);it.Valid();it.Next())
//        LOG::cout<<it.Key()<<" ->  "<<it.Data()->S<<std::endl;

}
//#####################################################################
// Function Set_Expression
//#####################################################################
template<class T> void DYNAMIC_IMPLICIT_SURFACE<T>::
Set_Expression(const std::string& expr)
{
    const char* str=expr.c_str();
    while(Parse(str)){}
}
//#####################################################################
// Function Principal_Curvatures
//#####################################################################
template<class T> VECTOR<T,2> DYNAMIC_IMPLICIT_SURFACE<T>::
Principal_Curvatures(const TV& X) const
{
    Eval(X);
    const TV& D=eval_list.Last()->D;
    const MATRIX<T,3>& H=eval_list.Last()->H;
    MATRIX<T,3> P=(T)1-MATRIX<T,3>::Outer_Product(H*D,D);
    MATRIX<T,3> B=P*H*P;
    T t=B.Trace()/2,t2=(B*B).Trace()/2;
    T d=sqrt(std::max((T)0,t2-t*t));
//    LOG::cout<<"curve: "<<X<<"  ->  "<<VECTOR<T,2>(t-d,t+d)<<std::endl;
    return VECTOR<T,2>(t-d,t+d);
}
template class DYNAMIC_IMPLICIT_SURFACE<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class DYNAMIC_IMPLICIT_SURFACE<double>;
#endif
}
