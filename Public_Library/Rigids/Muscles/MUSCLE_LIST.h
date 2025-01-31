//#####################################################################
// Copyright 2005-2007, Kevin Der, Eftychios Sifakis, Jonathan Su, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MUSCLE_LIST
//#####################################################################
#ifndef __MUSCLE_LIST__
#define __MUSCLE_LIST__

#include <Core/Arrays/ARRAY.h>
#include <Core/Data_Structures/TRIPLE.h>
#include <Rigids/Muscles/MUSCLE.h>
#include <Rigids/Muscles/MUSCLE_FORCE_CURVE.h>
namespace PhysBAM{

template<class TV> struct ATTACHMENT_POINT;
template<class TV> class RIGID_BODY_COLLECTION;
class VIEWER_DIR;

template<class TV>
class MUSCLE_LIST
{
    typedef typename TV::SCALAR T;
public:
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection;
    MUSCLE_FORCE_CURVE<T> muscle_force_curve;
    ARRAY<MUSCLE<TV>*> muscles;
    ARRAY<T> muscle_activations;
    ARRAY<ARRAY<TRIPLE<int,ATTACHMENT_POINT<TV>*,ATTACHMENT_POINT<TV>*> >,int> muscle_attachments_on_rigid_body;

    MUSCLE_LIST(RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input)
        :rigid_body_collection(rigid_body_collection_input)
    {}

    MUSCLE_LIST(const MUSCLE_LIST&) = delete;
    void operator=(const MUSCLE_LIST&) = delete;

    void Clean_Memory()
    {muscles.Delete_Pointers_And_Clean_Memory();muscle_activations.Clean_Memory();}

    int Add_Muscle(MUSCLE<TV>* new_muscle)
    {muscles.Append(new_muscle);new_muscle->id=muscles.m;muscle_activations.Append(1);return muscles.m;}

    void Set_Muscle_Activation(const int muscle_id,const T activation)
    {muscle_activations(muscle_id)=activation;}

//#####################################################################
    void Initialize_Muscle_Attachments_On_Rigid_Body();
    void Read(const VIEWER_DIR& viewer_dir);
    void Write(const STREAM_TYPE stream_type,const VIEWER_DIR& viewer_dir) const;
//##################################################################### 
};
}
#endif
