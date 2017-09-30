//#####################################################################
// Copyright 2017, Ounan Ding, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Log/LOG.h>
#include <Core/Math_Tools/RANGE_ITERATOR.h>
namespace PhysBAM{
int range_lookup[2][2]=
{
    // indices, side_adj
    {0x333000,0x888888}, // interior
    {0x312001,0x878888}  // ghost
};
//#####################################################################
// Constructor
//#####################################################################
template<int d> RANGE_ITERATOR<d>::
RANGE_ITERATOR(const RANGE<TV_INT>& outer,const RANGE<TV_INT>& inner,
    RI flags,int side_input)
    :vecs{outer.min_corner,inner.min_corner,inner.max_corner-1,outer.max_corner-1}
{
    Initialize(flags,side_input);
}
//#####################################################################
// Constructor
//#####################################################################
template<int d> RANGE_ITERATOR<d>::
RANGE_ITERATOR(const RANGE<TV_INT>& range,int outer_ghost,int inner_ghost,
    RI flags,int side_input)
    :vecs{range.min_corner-outer_ghost,range.min_corner-inner_ghost,
        range.max_corner+(inner_ghost-1),range.max_corner+(outer_ghost-1)}
{
    Initialize(flags,side_input);
}
//#####################################################################
// Constructor
//#####################################################################
template<int d> RANGE_ITERATOR<d>::
RANGE_ITERATOR(const TV_INT& counts,int outer_ghost,int inner_ghost,
    RI flags,int side_input)
    :vecs{TV_INT()-outer_ghost,TV_INT()-inner_ghost,
        counts+(inner_ghost-1),counts+(outer_ghost-1)}
{
    Initialize(flags,side_input);
}
//#####################################################################
// Constructor
//#####################################################################
template<int d> RANGE_ITERATOR<d>::
RANGE_ITERATOR(const RANGE<TV_INT>& range,RI flags)
    :vecs{range.min_corner,range.min_corner,range.max_corner-1,range.max_corner-1}
{
    Initialize(flags,-1);
}
//#####################################################################
// Constructor
//#####################################################################
template<int d> RANGE_ITERATOR<d>::
RANGE_ITERATOR(const TV_INT& counts,int outer_ghost,RI flags)
    :vecs{TV_INT()-outer_ghost,TV_INT()-outer_ghost,
        counts+(outer_ghost-1),counts+(outer_ghost-1)}
{
    Initialize(flags,-1);
}
//#####################################################################
// Function Set_Range
//#####################################################################
template<int d> void RANGE_ITERATOR<d>::
Set_Range(const RANGE<TV_INT>& outer,const RANGE<TV_INT>& inner)
{
    vecs[0]=outer.min_corner;
    vecs[1]=inner.min_corner;
    vecs[2]=inner.max_corner-1;
    vecs[3]=outer.max_corner-1;
}
//#####################################################################
// Function Set_Range
//#####################################################################
template<int d> void RANGE_ITERATOR<d>::
Set_Range(const RANGE<TV_INT>& range,int outer_ghost,int inner_ghost)
{
    vecs[0]=range.min_corner-outer_ghost;
    vecs[1]=range.min_corner-inner_ghost;
    vecs[2]=range.max_corner+(inner_ghost-1);
    vecs[3]=range.max_corner+(outer_ghost-1);
}
//#####################################################################
// Function Set_Range
//#####################################################################
template<int d> void RANGE_ITERATOR<d>::
Set_Range(const TV_INT& counts,int outer_ghost,int inner_ghost)
{
    vecs[0]=TV_INT()-outer_ghost;
    vecs[1]=TV_INT()-inner_ghost;
    vecs[2]=counts+(inner_ghost-1);
    vecs[3]=counts+(outer_ghost-1);
}
//#####################################################################
// Function Initialize
//#####################################################################
template<int d> void RANGE_ITERATOR<d>::
Initialize(RI flags,int side_input)
{
    Encode(flags,side_input);
    Reset(flags);
    if(any(flags&RI::end)) index(d-1)++;
}
//#####################################################################
// Function Reset
//#####################################################################
template<int d> void RANGE_ITERATOR<d>::
Reset(RI flags)
{
    if(!(flags&(RI::reverse|RI::end))){
        side=-1;
        Next_Side();}
    else{
        side=2*d;
        Prev_Side();}
}
//#####################################################################
// Function Encode
//#####################################################################
template<int d> void RANGE_ITERATOR<d>::
Encode(RI flags,int side_input)
{
    bool one_side=!(flags&RI::side_mask) && side_input>=0;
    if(any(flags&RI::side_mask)) side_mask=side_input&((1<<2*d)-1);
    else side_mask=side_input<0?(1<<2*d)-1:(1<<side_input);
    flags=flags&~(RI::reverse|RI::end);
    if(!any(flags&RI::ghost)) side_mask=1;

    // Explicitly check for flag combinations that do not make sense.
    if(!any(flags&RI::ghost)){
        PHYSBAM_ASSERT(!(flags&~(RI::ghost)));}
    if(any(flags&RI::duplicate_corners)){
        PHYSBAM_ASSERT(!(flags&~(RI::ghost|RI::duplicate_corners|RI::side_mask)));}

    int offsets[2];
    int mode=(int)(flags&RI::ghost);
    for(int i=0;i<2;i++) offsets[i]=range_lookup[mode][i];

    if(one_side && !(flags&RI::partial_single_side))
        for(int i=0;i<2;i++)
            offsets[i]=((offsets[i]&0xF00F00)>>8)|(offsets[i]&0xFF0FF0);

    if(any(flags&RI::delay_corners))
        for(int i=0;i<2;i++)
            offsets[i]=((offsets[i]&0x00F00F)<<8)|((offsets[i]&0xF00F00)>>8)|(offsets[i]&0x0F00F0);

    if(any(flags&RI::duplicate_corners))
        for(int i=0;i<2;i++)
            offsets[i]=((offsets[i]&0xF00F00)>>8)|(offsets[i]&0xFF0FF0);

    indices=offsets[0];
    side_adj=offsets[1];
}
//#####################################################################
// Function Next_Helper
//#####################################################################
template<int d> void RANGE_ITERATOR<d>::
Next_Helper()
{
    // next row of range
    index(d-1)=current[0](d-1);
    for(int i=d-2;i>=0;i--){
        if(++index(i)<current[1](i)) return;
        index(i)=current[0](i);}
    Next_Side();
}
//#####################################################################
// Function Prev_Helper
//#####################################################################
template<int d> void RANGE_ITERATOR<d>::
Prev_Helper()
{
    // next row of range
    index(d-1)=current[1](d-1)-1;
    for(int i=d-2;i>=0;i--){
        if(--index(i)>=current[0](i)) return;
        index(i)=current[1](i)-1;}
    Prev_Side();
}
//#####################################################################
// Function Next_Side
//#####################################################################
template<int d> void RANGE_ITERATOR<d>::
Next_Side()
{
    // finished current block; next block
    side++;
    if(!(side_mask>>side)){
        index(d-1)=current[1](d-1); // Make invalid
        return;} // Done
    while(!(side_mask&(1<<side))) side++;
    Fill_Current();
    index=current[0];
    if(!current[0].All_Less(current[1])) Next_Side();
}
//#####################################################################
// Function Prev_Side
//#####################################################################
template<int d> void RANGE_ITERATOR<d>::
Prev_Side()
{
    // finished current block; next block
    side--;
    if(side<0 || !(side_mask&((2<<side)-1))){
        index(d-1)=current[0](d-1)-1; // Make invalid
        return;} // Done
    while(!(side_mask&(1<<side))) side--;
    Fill_Current();
    index=current[1]-1;
    if(!current[0].All_Less(current[1])) Prev_Side();
}
//#####################################################################
// Function Fill_Current
//#####################################################################
template<int d> void RANGE_ITERATOR<d>::
Fill_Current()
{
    // Assemble next block from rules.
    int side_axis=side/2,s1=side&1;
    for(int i=0;i<d;i++){
        int p=(i>=side_axis)+(i>side_axis);
        for(int s=0;s<2;s++){
            int sh=(3*s+p)*4,k=s;
            int j=(indices>>sh)&3;
            int o=((side_adj>>sh)&15)-8;
            if(s1 && i==side_axis){j^=3;o=-o;k^=1;}
            current[k](i)=vecs[j](i)+o;}}
    current[1]+=1;
}
//#####################################################################
template class RANGE_ITERATOR<0>;
template class RANGE_ITERATOR<1>;
template class RANGE_ITERATOR<2>;
template class RANGE_ITERATOR<3>;
}
