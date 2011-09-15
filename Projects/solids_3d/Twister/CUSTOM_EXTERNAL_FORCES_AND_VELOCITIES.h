//#####################################################################
// Copyright 2002, Ronald Fedkiw, Robert Bridson.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CUSTOM_EXTERNAL_FORCES_AND_VELOCITIES
//##################################################################### 
//
//#####################################################################
// Fedkiw - August 9, 2002
// Bridson - May 21, 2003
//#####################################################################
#ifndef __CUSTOM_EXTERNAL_FORCES_AND_VELOCITIES__
#define __CUSTOM_EXTERNAL_FORCES_AND_VELOCITIES__

namespace PhysBAM{

class CUSTOM_EXTERNAL_FORCES_AND_VELOCITIES:public EXTERNAL_FORCES_AND_VELOCITIES
{
public:
    int number_side_panels;
    double aspect_ratio;

    CUSTOM_EXTERNAL_FORCES_AND_VELOCITIES(const int number_side_panels_input,
                                                                                 const double aspect_ratio_input)
                                                                                 :number_side_panels(number_side_panels_input),aspect_ratio(aspect_ratio_input)
    {}

    virtual ~CUSTOM_EXTERNAL_FORCES_AND_VELOCITIES()
    {}

//#####################################################################
// Function Set_External_Velocities
//#####################################################################
void Set_External_Velocities(ARRAY<VECTOR_3D,VECTOR<int,1> >& V,const double time)
{
    // shape of cloth mesh
    int m=(int)(aspect_ratio*number_side_panels)+1,n=number_side_panels+1;
    
    int i,j;
    i=1;j=1;V(i+m*(j-1)).x=0;V(i+m*(j-1)).y=0;
    i=1;j=n;V(i+m*(j-1)).x=0;V(i+m*(j-1)).y=0;
    i=m;j=1;V(i+m*(j-1)).x=0;V(i+m*(j-1)).y=0;
    i=m;j=n;V(i+m*(j-1)).x=0;V(i+m*(j-1)).y=0;
}
//#####################################################################
// Function Zero_Out_Enslaved_Velocity_Nodes
//#####################################################################
void Zero_Out_Enslaved_Velocity_Nodes(ARRAY<VECTOR_3D,VECTOR<int,1> >& V,const double time)
{
    // shape of cloth mesh
    int m=(int)(aspect_ratio*number_side_panels)+1,n=number_side_panels+1;
    
    int i,j;
    i=1;j=1;V(i+m*(j-1)).x=0;V(i+m*(j-1)).y=0;
    i=1;j=n;V(i+m*(j-1)).x=0;V(i+m*(j-1)).y=0;
    i=m;j=1;V(i+m*(j-1)).x=0;V(i+m*(j-1)).y=0;
    i=m;j=n;V(i+m*(j-1)).x=0;V(i+m*(j-1)).y=0;
}
//#####################################################################
};
}
#endif




