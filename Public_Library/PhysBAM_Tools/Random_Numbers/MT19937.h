//#####################################################################
// Copyright 2002-2009, Wen Zheng.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MT19937
//#####################################################################
// Mersenne Twister pseudorandom number generator, MT19937 version.
//#####################################################################
#ifndef __MT19937__
#define __MT19937__

#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE_FUNCTIONS.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <iostream>
namespace PhysBAM{

template<class T>
class MT19937
{
public:
    typedef int HAS_UNTYPED_READ_WRITE;
    unsigned int mt[624];
    int index;
    static const int n=624,m=397;
    static const unsigned int a=0x9908b0df,UPPER_MASK=0x80000000,LOWER_MASK=0x7fffffff;

    MT19937();
    explicit MT19937(const unsigned int value);
    virtual ~MT19937();
    void Set_Seed(const unsigned int value=5489);
    T operator()(); // in [0,1)

    template<class RW> void Read(std::istream& input)
    {Read_Binary<RW>(input,index);Read_Binary_Array<RW>(input,mt,n);}

    template<class RW> void Write(std::ostream& output) const
    {Write_Binary<RW>(output,index);Write_Binary_Array<RW>(output,mt,n);}
//#####################################################################
};
template<class T> inline std::istream& 
operator>>(std::istream& input,MT19937<T>& g)
{input>>g.index;for(int i=0;i<g.n;i++) input>>g.mt[i];return input;}

template<class T> inline std::ostream&
operator<<(std::ostream& output,const MT19937<T>& g)
{output<<g.index;for(int i=0;i<g.n;i++) output<<" "<<g.mt[i];return output;}
}
#endif
