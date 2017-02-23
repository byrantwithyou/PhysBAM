//#####################################################################
// Copyright 2002-2009, Wen Zheng.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MT19937
//#####################################################################
#include <Core/Random_Numbers/MT19937.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class T> MT19937<T>::
MT19937(const unsigned int value)
{
    Set_Seed(value);
}
//#####################################################################
// Destructor
//#####################################################################
template<class T> MT19937<T>::
~MT19937()
{
}
//#####################################################################
// Function Set_Seed
//#####################################################################
template<class T> void MT19937<T>::
Set_Seed(const unsigned int value)
{
    mt[0]=value&0xffffffffUL;
    for(index=1;index<n;index++){
        mt[index]=(1812433253UL*(mt[index-1]^(mt[index-1]>>30))+index); 
        mt[index]&=0xffffffffUL;}
}
//#####################################################################
// Operator ()
//#####################################################################
template<class T> T MT19937<T>::
operator()()
{
    unsigned int y;
    static unsigned int mag01[2]={0x0UL,a};

    if(index>=n){
        int kk;
        for(kk=0;kk<n-m;kk++){
            y=(mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk]=mt[kk+m]^(y>>1)^mag01[y&0x1UL];}
        for(;kk<n-1;kk++){
            y=(mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk]=mt[kk+(m-n)]^(y>>1)^mag01[y&0x1UL];}
        y=(mt[n-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        mt[n-1]=mt[m-1]^(y>>1)^mag01[y&0x1UL];
        index=0;}
  
    y=mt[index++];
    y^=(y>>11);
    y^=(y<<7)&0x9d2c5680UL;
    y^=(y<<15)&0xefc60000UL;
    y^=(y>>18);
    return T(y*(1.0/4294967296.0));
}
//#####################################################################
template class MT19937<double>;
template class MT19937<float>;
}
