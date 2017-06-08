//#####################################################################
// Copyright 2017, Ounan Ding, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Log/LOG.h>
#include <Grid_Tools/Grids/FACE_RANGE_ITERATOR.h>
namespace PhysBAM{
int face_range_lookup[2][5]=
{
    // indices, axis_adj, side_adj, skip_outer, skip_inner
    {0x312001,0x887889,0x878888,0x788998,0x879887}, // NO FLAGS
    {0x333000,0x888888,0x888888,0x777999,0x888888}  // FLAG_INTERIOR
};
//#####################################################################
// Constructor
//#####################################################################
template<int d> FACE_RANGE_ITERATOR<d>::
FACE_RANGE_ITERATOR(const RANGE<TV_INT>& outer,const RANGE<TV_INT>& inner,
    RF flags,int axis,int side_input)
    :vecs{outer.min_corner,inner.min_corner,inner.max_corner-1,outer.max_corner-1}
{
    Initialize(axis,side_input,flags);
}
//#####################################################################
// Constructor
//#####################################################################
template<int d> FACE_RANGE_ITERATOR<d>::
FACE_RANGE_ITERATOR(const RANGE<TV_INT>& range,int outer_ghost,int inner_ghost,
    RF flags,int axis,int side_input)
    :vecs{range.min_corner-outer_ghost,range.min_corner-inner_ghost,
        range.max_corner+(inner_ghost-1),range.max_corner+(outer_ghost-1)}
{
    Initialize(axis,side_input,flags);
}
//#####################################################################
// Constructor
//#####################################################################
template<int d> FACE_RANGE_ITERATOR<d>::
FACE_RANGE_ITERATOR(const TV_INT& counts,int outer_ghost,int inner_ghost,
    RF flags,int axis,int side_input)
    :vecs{TV_INT()-outer_ghost,TV_INT()-inner_ghost,
        counts+(inner_ghost-1),counts+(outer_ghost-1)}
{
    Initialize(axis,side_input,flags);
}
//#####################################################################
// Constructor
//#####################################################################
template<int d> FACE_RANGE_ITERATOR<d>::
FACE_RANGE_ITERATOR(const RANGE<TV_INT>& range,
    RF flags,int axis,int side_input) // implict RF::interior
    :vecs{range.min_corner,range.min_corner,range.max_corner-1,range.max_corner-1}
{
    Initialize(axis,side_input,flags|RF::interior);
}
//#####################################################################
// Constructor
//#####################################################################
template<int d> FACE_RANGE_ITERATOR<d>::
FACE_RANGE_ITERATOR(const TV_INT& counts,int outer_ghost,
    RF flags,int axis,int side_input) // implict RF::interior
    :vecs{TV_INT()-outer_ghost,TV_INT()-outer_ghost,
        counts+(outer_ghost-1),counts+(outer_ghost-1)}
{
    Initialize(axis,side_input,flags|RF::interior);
}
//#####################################################################
// Function Encode
//#####################################################################
template<int d> void FACE_RANGE_ITERATOR<d>::
Initialize(int axis,int side_input,RF flags)
{
    Encode(axis,side_input,flags);
    Reset(flags);
    if(any(flags&RF::end)) face.index(d-1)++;
}
//#####################################################################
// Function Reset_Axis_Side
//#####################################################################
template<int d> void FACE_RANGE_ITERATOR<d>::
Reset(RF flags)
{
    if(!(flags&(RF::reverse|RF::end))){
        side=2*d;
        face.axis=-1;
        Next_Axis_Side();}
    else{
        side=-1;
        face.axis=d;
        Prev_Axis_Side();}
}
//#####################################################################
// Function Encode
//#####################################################################
template<int d> void FACE_RANGE_ITERATOR<d>::
Encode(int axis,int side_input,RF flags)
{
    bool one_side=!(flags&RF::side_mask) && side_input>=0;
    if(any(flags&RF::axis_mask)) axis_mask=axis;
    else axis_mask=axis<0?(1<<d)-1:(1<<axis);
    if(any(flags&RF::side_mask)) side_mask=side_input;
    else side_mask=side_input<0?(1<<2*d)-1:(1<<side_input);
    flags=flags&~(RF::reverse|RF::end|RF::axis_mask);
    if(any(flags&RF::interior)) side_mask=1;

    // Explicitly check for flag combinations that do not make sense.
    if(any(flags&RF::interior)){
        PHYSBAM_ASSERT(!(flags&~(RF::interior|RF::skip_outer)));}
    if(any(flags&RF::duplicate_corners)){
        PHYSBAM_ASSERT(!(flags&~(RF::duplicate_corners|RF::skip_outer|RF::skip_inner|RF::side_mask)));}

    int offsets[3];
    int mode=(int)(flags&RF::interior);
    for(int i=0;i<3;i++) offsets[i]=face_range_lookup[mode][i];
    if(any(flags&RF::skip_outer)) offsets[1]+=face_range_lookup[mode][3]-0x888888;
    if(any(flags&RF::skip_inner)) offsets[1]+=face_range_lookup[mode][4]-0x888888;

    if(one_side && !(flags&RF::partial_single_side))
        for(int i=0;i<3;i++)
            offsets[i]=((offsets[i]&0xF00F00)>>8)|(offsets[i]&0xFF0FF0);

    if(any(flags&RF::delay_corners))
        for(int i=0;i<3;i++)
            offsets[i]=((offsets[i]&0x00F00F)<<8)|((offsets[i]&0xF00F00)>>8)|(offsets[i]&0x0F00F0);

    if(any(flags&RF::duplicate_corners))
        for(int i=0;i<3;i++)
            offsets[i]=((offsets[i]&0xF00F00)>>8)|(offsets[i]&0xFF0FF0);

    indices=offsets[0];
    axis_adj=offsets[1];
    side_adj=offsets[2];
}
//#####################################################################
// Function Next_Helper
//#####################################################################
template<int d> void FACE_RANGE_ITERATOR<d>::
Next_Helper()
{
    // next row of range
    face.index(d-1)=current[0](d-1);
    for(int i=d-2;i>=0;i--){
        if(++face.index(i)<current[1](i)) return;
        face.index(i)=current[0](i);}
    Next_Axis_Side();
}
//#####################################################################
// Function Next_Helper
//#####################################################################
template<int d> void FACE_RANGE_ITERATOR<d>::
Prev_Helper()
{
    // next row of range
    face.index(d-1)=current[1](d-1)-1;
    for(int i=d-2;i>=0;i--){
        if(--face.index(i)>=current[0](i)) return;
        face.index(i)=current[1](i)-1;}
    Prev_Axis_Side();
}
//#####################################################################
// Function Next_Helper
//#####################################################################
template<int d> void FACE_RANGE_ITERATOR<d>::
Next_Axis_Side()
{
    // finished current block; next block
    side++;
    if(!(side_mask>>side)){
        side=0;
        face.axis++;
        if(!(axis_mask>>face.axis)){
            face.index(d-1)=current[1](d-1); // Make invalid
            return;} // Done
        while(!(axis_mask&(1<<face.axis))) face.axis++;}
    while(!(side_mask&(1<<side))) side++;
    Fill_Current();
    face.index=current[0];
    if(!current[0].All_Less(current[1])) Next_Axis_Side();
}
//#####################################################################
// Function Next_Helper
//#####################################################################
template<int d> void FACE_RANGE_ITERATOR<d>::
Prev_Axis_Side()
{
    // finished current block; next block
    side--;
    if(side<0 || !(side_mask&((2<<side)-1))){
        side=2*d-1;
        face.axis--;
        if(face.axis<0 || !(axis_mask&((2<<face.axis)-1))){
            face.index(d-1)=current[0](d-1)-1; // Make invalid
            return;} // Done
        while(!(axis_mask&(1<<face.axis))) face.axis--;}
    while(!(side_mask&(1<<side))) side--;
    Fill_Current();
    face.index=current[1]-1;
    if(!current[0].All_Less(current[1])) Prev_Axis_Side();
}
//#####################################################################
// Function Fill_Current
//#####################################################################
template<int d> void FACE_RANGE_ITERATOR<d>::
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
    int p=(face.axis>=side_axis)+(face.axis>side_axis);
    for(int s=0;s<2;s++){
        int sh=(3*s+p)*4,k=s;
        int o=((axis_adj>>sh)&15)-8;
        if(s1 && face.axis==side_axis){o=-o;k^=1;}
        current[k](face.axis)+=o;}
    current[1](face.axis)++;
    current[1]+=1;
    PHYSBAM_ASSERT(current[0].All_Less_Equal(current[1]));
}
//#####################################################################
template class FACE_RANGE_ITERATOR<1>;
template class FACE_RANGE_ITERATOR<2>;
template class FACE_RANGE_ITERATOR<3>;
}
