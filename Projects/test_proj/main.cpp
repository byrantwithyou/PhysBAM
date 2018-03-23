//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include <Core/Log/LOG.h>
#include <Core/Vectors/VECTOR.h>

using namespace PhysBAM;

typedef float RW;
typedef double T;
typedef VECTOR<T,3> TV;
typedef VECTOR<int,TV::m> TV_INT;

template<class... Args>
void test_chr(char* buff,Args... args)
{
    printf("%s %s",__FUNCTION__,buff);
    int a=LOG::printf(buff,17,args...,'q',19,1.3);
    int b=printf(buff,17,args...,'q',19,1.3);
    printf("# %i %i\n",a,b);
}

template<class... Args>
void test_cpt(char* buff,Args... args)
{
    printf("%s %s",__FUNCTION__,buff);
    int r=99;
    int a=LOG::printf(buff,17,args...,&r,19,1.3);
    int b=printf(buff,17,args...,&r,19,1.3);
    printf("# %i %i\n",a,b);
}


template<class... Args>
void test_flt(char* buff,Args... args)
{
    printf("%s %s",__FUNCTION__,buff);
    int a=LOG::printf(buff,17,args...,7.1,19,1.3);
    int b=printf(buff,17,args...,7.1,19,1.3);
    printf("# %i %i\n",a,b);
}


template<class... Args>
void test_int(char* buff,Args... args)
{
    printf("%s %s",__FUNCTION__,buff);
    int a=LOG::printf(buff,17,args...,13,19,1.3);
    int b=printf(buff,17,args...,13,19,1.3);
    printf("# %i %i\n",a,b);
}


template<class... Args>
void test_ptr(char* buff,Args... args)
{
    printf("%s %s",__FUNCTION__,buff);
    int a=LOG::printf(buff,17,args...,(void*)buff,19,1.3);
    int b=printf(buff,17,args...,(void*)buff,19,1.3);
    printf("# %i %i\n",a,b);
}


template<class... Args>
void test_str(char* buff,Args... args)
{
    printf("%s %s",__FUNCTION__,buff);
    int a=LOG::printf(buff,17,args...,"xyz",19,1.3);
    int b=printf(buff,17,args...,"xyz",19,1.3);
    printf("# %i %i\n",a,b);
}

template<class... Args>
void gen_type(char* buff,int n,Args... args)
{
    strcpy(buff+n+1,"B%iC%gD\n");
    buff[n]='d';
    test_int(buff,args...);
    buff[n]='i';
    test_int(buff,args...);
    buff[n]='u';
    test_int(buff,args...);
    buff[n]='s';
    test_str(buff,args...);
    buff[n]='c';
    test_chr(buff,args...);
    buff[n]='p';
    test_ptr(buff,args...);
    buff[n]='F';
    test_flt(buff,args...);
    buff[n]='f';
    test_flt(buff,args...);
    buff[n]='E';
    test_flt(buff,args...);
    buff[n]='e';
    test_flt(buff,args...);
    buff[n]='G';
    test_flt(buff,args...);
    buff[n]='g';
    test_flt(buff,args...);
    buff[n]='X';
    test_int(buff,args...);
    buff[n]='x';
    test_int(buff,args...);
    buff[n]='o';
    test_int(buff,args...);
    // buff[n]='A';
    // test_flt(buff,args...);
    // buff[n]='a';
    // test_flt(buff,args...);
    buff[n]='n';
    test_cpt(buff,args...);
}

template<class... Args>
void gen_prec(char* buff,int n,Args... args)
{
    gen_type(buff,n,args...);
    buff[n++]='.';
    buff[n]='*';
    gen_type(buff,n+1,args...,8);
    buff[n]='9';
    gen_type(buff,n+1,args...);
}

void gen_width(char* buff,int n)
{
    gen_prec(buff,n);
    buff[n]='*';
    gen_prec(buff,n+1,4);
    buff[n]='5';
    gen_prec(buff,n+1);
}


void gen_format(char* buff,int n)
{
    buff[n++]='%';
    for(int i=0;i<32;i++)
    {
        int k=n;
        if(i&1) buff[k++]='+';
        if(i&2) buff[k++]='-';
        if(i&4) buff[k++]=' ';
        if(i&8) buff[k++]='0';
        if(i&16) buff[k++]='#';
        gen_width(buff,k);
    }
}

int main(int argc, char* argv[])
{
    printf("%5.3i\n",2);
    printf("%3.5i\n",2);


    char buff[1000]="%iA";
    gen_format(buff,3);



    
    return 0;
}

