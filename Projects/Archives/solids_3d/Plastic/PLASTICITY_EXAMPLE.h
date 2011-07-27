//#####################################################################
// Copyright 2003, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PLASTIC_SPHERE_EXAMPLE
//#####################################################################
#ifndef __PLASTICITY_EXAMPLE__
#define __PLASTICITY_EXAMPLE__

#include "TET_SIM_EXAMPLE.h"
#include <Forces_And_Torques/DIAGONALIZED_PLASTICITY_CONTROL_3D.h>
#include <Forces_And_Torques/DIAGONALIZED_SIMPLE_PLASTICITY_3D.h>
namespace PhysBAM{

template<class T,class RW>
class PLASTICITY_EXAMPLE:public TET_SIM_EXAMPLE<T,RW>
{
public:
    bool use_plasticity,use_control;
    bool preserve_volume,show_goal;
    DIAGONALIZED_PLASTICITY_MODEL_3D<T>* plasticity_model;
    T yield_ratio,plastic_clamp_ratio;
    ARRAY<VECTOR_3D<T> > plastic_goal;
    char input_file[256];

    PLASTICITY_EXAMPLE()
        :use_plasticity(false),use_control(false),preserve_volume(false),plasticity_model(0),yield_ratio(2),plastic_clamp_ratio(2)
    {}

    ~PLASTICITY_EXAMPLE()
    {}
    
    VECTOR_3D<T> Linear_Map(VECTOR_3D<T>& p)
    {SYMMETRIC_MATRIX<T,3> G=SYMMETRIC_MATRIX<T,3>::Identity_Matrix();
    G-=(T).75*SYMMETRIC_MATRIX<T,3>::Outer_Product(VECTOR_3D<T>(1,(T)1.5,0).Normalized());
    return G*p;}
    
    VECTOR_3D<T> Wave_Map(VECTOR_3D<T>& p)
    {VECTOR_3D<T> a1=(T).006*VECTOR_3D<T>(1,1,1),a2=(T).008*VECTOR_3D<T>(2,3,-4);
    VECTOR_3D<T> d1=(T)3*VECTOR_3D<T>(2,3,4),d2=(T)5*VECTOR_3D<T>(-2,3,1);
    return p+sin(VECTOR_3D<T>::Dot_Product(d1,p))*a1+sin(VECTOR_3D<T>::Dot_Product(d2,p))*a2;}

    VECTOR_3D<T> Bump_Map(VECTOR_3D<T>& p)
    {VECTOR_3D<T> d=p.Normalized();
    T scale=1+(T).1*sin(6*d.x)*sin(9*d.y)*sin(15*d.z);
    return scale*p;}

    VECTOR_3D<T> Pancake_Map(VECTOR_3D<T>& p,const T goal_thickness,const T threshold,const MATRIX<T,3>& Q)
    {VECTOR_3D<T> pr=Q*p;
    T thickness=sqrt(1-sqr(pr.y)-sqr(pr.z));
    T t=thickness/threshold;
    T new_thickness=goal_thickness*(t>1?1:t*(2-t));
    T x_scale=new_thickness/thickness;
    return Q.Transposed()*VECTOR_3D<T>(x_scale*pr.x,pr.y,pr.z);}
    
//#####################################################################
// Function Initialize_Diagonalized_Finite_Volume_Model
//#####################################################################
void Initialize_Diagonalized_Finite_Volume_Model(DEFORMABLE_OBJECT<T,VECTOR_3D<T> >& deformable_object,TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume)
{
    TET_SIM_EXAMPLE<T,RW>::Initialize_Diagonalized_Finite_Volume_Model(deformable_object,tetrahedralized_volume);
    if(use_control){
        DIAGONALIZED_PLASTICITY_CONTROL_3D<T>* control=new DIAGONALIZED_PLASTICITY_CONTROL_3D<T>(*strain,plastic_goal,yield_ratio);
        control->Print_Extreme_Deformations();plasticity_model=control;}
    else if(use_plasticity)plasticity_model=new DIAGONALIZED_SIMPLE_PLASTICITY_3D<T>(strain->Dm_inverse.m,yield_ratio,plastic_clamp_ratio);
    diagonalized_fvm->plasticity_model=plasticity_model;
    plastic_goal.Clean_Memory();
}
//#####################################################################
// Function Read_Data_Files
//#####################################################################
void Read_Data_Files(TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume,T& time,const int frame)
{
    TET_SIM_EXAMPLE<T,RW>::Read_Data_Files(tetrahedralized_volume,time,frame);
    if(plasticity_model)plasticity_model->template Read_State<RW>(output_directory,frame);
}
//#####################################################################
// Function Write_Data_Files
//#####################################################################
void Write_Data_Files(const TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume,const ARRAY<RIGID_BODY<TV>*>& rigid_bodies,const T time,const int frame)
{
    TET_SIM_EXAMPLE<T,RW>::Write_Data_Files(tetrahedralized_volume,rigid_bodies,time,frame);
    if(plasticity_model)plasticity_model->template Write_State<RW>(output_directory,frame);
}
//#####################################################################
};
}
#endif
